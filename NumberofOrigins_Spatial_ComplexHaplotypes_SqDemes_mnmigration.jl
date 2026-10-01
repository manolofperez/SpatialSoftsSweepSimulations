using PoissonRandom
using Distributions
using LinearAlgebra
using StatsBase
using MAT
using Luxor
using Plots
using Colors, ColorSchemes

include("GaussPoissonHybrid_mnrnd.jl")#consider removing this as it is typically slower (somehow) than Julia's multinomial random number generator



@doc """
Main simulation function to simulate the infinite alleles process of many origins of a selective allele, which includes genetic drift, mutation, selection and migration
for a population of size N. Local migration is implemented assuming effectively continous space diffusion, by having a fixed migration rate p to adjacent demes, 
but with deme sizes dependent on the diffusion constant. 

The simulation then samples the true population at specific locations and at specific times corresponding to how populations were sampled in field and produces 
a single array containing sampled alleles

Input:
- N: population size of chromosomes (if simulating diploid organism with additive fitness effects divide by 2),
- Area: Total area (this is ignored if shape outline specified in lon and lat aren't empty),
- theta: population-scaled mutation rate (2Nμ),
- s: selection coefficient,
- D: diffusion constant [km^2 per generation]
- Ratio: aspect ratio (width / height) (this is ignored if shape outline specified in lon and lat aren't empty)
- Gaussian: 1 (using Gaussian approximation) or 0 (multinomial distribution)
- maxt: number of generations to run simulations
- xmut_target (if this is not empty then simulations end when mutant frequency reaches this value, instead of at maxt)
- lon and lat contain arrays defining the species range in longitude and latitude
- Nsamp: number of chromosomes/haploid individuals sampled at each spatial location 
- samplingtime:
- Nsampt: number of chromosomes/haploid individuals sampled at each space-time location
- uniquespatialindx: index of length Nsampt into Nsamp array indicating the unique spatial locations
- spacetimesample_deme_number: the number of the nearest deme for each spatial sampling location
- uniquesamplingtime:

Output:
- Xt: cell array of length = number of demes containing haplotype number and frequency array in that deme,
      stored as a complex number h + i*f
- ncomplexhap: cell array of length number of spatial sampling sites containing haplotype number 
    and *number* of chromosomes/haploids sampled in that deme, stored as complex number h +i*n
- nx: Number of demes in x direction
- Ademe: Area of each deme
- NDeme: Number of haploid individuals/chromosomes per deme
- t: final generation time at end of simulation
- truefreq: true frequency of mutant allele

"""
function NumberOfOrigins_Spatial(N::Real, Area::Real, theta::Real, s::Real, D::Real, hdom::Real, Ratio::Real, Gaussian::Real, maxt::Real, xmut_target, lon,lat, Nsamp,samplingtime,Nsampt,uniquespatialindx,spacetimesample_deme_number,uniquesamplingtime)
    

    alpha = 2*N*s
    p=0.05 # percentage that migrate out of deme per generation — keep fixed at p=0.05 unless you have good reason!
    println("")
    println("")
    println("N = ",round(N/1e6;sigdigits=3)," million")
    println("s = ",round(s;sigdigits = 3))
    println("D = ",round(D;sigdigits = 3))

    
    # println("Area = ",Area)
    println("2Nμ = ",round(theta;sigdigits = 2))
    println("2Ns = ",round(alpha;sigdigits = 2))

    local μ = theta / 2 / N # mutation
    println("μ = ",μ)

    println("Fixed probability of migration p = ",p)
    
    println("Hybrid Gauss-Poisson approx = ",Gaussian)
     
    println("max time (time in generations from introduction of insecticides to today) = ", maxt)

    if !isempty(xmut_target)
        # xmut_target = 0.692 #set this to [] for normal simulations that end at t=maxt
        println("Target  mutant sample frequency = ", round(xmut_target;sigdigits=3))
    end
    #in this case sampling is done exactly at the times and locations sampled as specified in samplingtime and spacetimesample_deme_number,
    #otherwise if xmut_target is not empty simulations end when threshold sample frequency of mutants is reached and sampling effectively 
    #ignores time and for each deme the number of chromosomes sampled every generation is the number specified in Nsampt, where we collect 
    #samples over multiple sample sites in the same deme and over the corresponding multiple collection times for those sites


    maxt = maxt + 1 #time index in loop below start at t=1, so need to add one to maxt

    
    #Deme geographical parameters
    Ademe = 4*D/p #square demes, where we assume p*w^2 = 4D is the MSD in one generation and w is the deme width
    wdeme = sqrt(Ademe)

    println("Area of each deme  = ",round(Ademe;sigdigits=3),"km^2")
    println("Width of each square deme  = ",round(wdeme;sigdigits=3),"km")

    truefreq = [0.0] #mutant frequency
 
    println("uniquesamplingtime, =", uniquesamplingtime)

    #convert shape file inputs of lon and lat in degrees to x and y in kms
    R=6378.137#km radius at the equator — defined here (not inside the if-block below) so it's always available, including when lon is empty

    if !isempty(lon)

  
        #calculate x,y positions of shape outline in lon and lat

        ND, nx, ny, x, y, minx, miny = CalculateDemes(Ademe, Area, Ratio,lon,lat,R) #x and y have their minimum subtracted

      

        
    else


        ND, nx, ny = CalculateDemes(Ademe, Area, Ratio, lon,lat,R)


    end


    
    #Create CartesianIndices over deme array of size ny x nx
    CartInd = CartesianIndices((1:ny,1:nx))

    mutant_counter = 0 #Counts the number of de novo mutations that occur in the simulation


    #Intial guess of max number of origins in each deme
    maxη=40

    KK = ceil(Int,maxη)
    K = [KK for ii in eachindex(CartInd)]


    if !isempty(lon)
        #Create a mask array for Africa/sub-Saharan Africa from shape outline input in x and y in kms
        polygonAfrica = Point.(vec(x),vec(y)) #create polygon object from x and y positions defining outline/shape of Africa

        xp = 0:wdeme:(nx-1)*wdeme #vector of deme positions/points in x-direction
        yp = 0:wdeme:(ny-1)*wdeme #vector of deme positions/points in y-direction
        demepoints = vec(Point.(xp',yp)) #calling Point in this way, where one argument is row vector and the other column — like an outer product — gives a *vector* of deme positions: size(demepoints) = nx*ny
        demepoints = vec(demepoints)
        # println("size(demepoints) = ", size(demepoints))
        DemesInAfrica = [isinside(pts,polygonAfrica; allowonedge=true) for pts in demepoints] #size(DemesInAfrica) = nx*ny — are these column order and which way?
        #Both demepoints and DemesInAfrica are column-ordered with increasing y-values, then increasing x values:
        #e.g. demepoints = (0,0) , (0,wdeme), (0,2*wdeme), ... (0, (nx-1)*wdeme), (wdeme,0), (wdeme, wdeme), (wdeme, 2*wdeme),...
        # => my arbirary integer grid points for demes will ba mapped onto space in this way

        #Calculate area from polyarea in Luxor and shape file of gambiae range
        Area = polyarea(polygonAfrica)


    else
        DemesInAfrica = fill(true,ND)

    end

    #Calculate nearest neighbours for each deme and update DemesInAfrica to remove any demes with zero nearest neighbours
    nnind, nnn, DemesInAfrica = CalcNearestNeighbours(nx,ny, DemesInAfrica)

    ND = length(findall(DemesInAfrica))#actual number of demes in sub-Saharan Africa
    NDeme = round(Int, N / ND) #Population of mosquitos in a Deme 
    println("Number of demes ND = ",ND)
    println("Population size per deme NDeme = ", NDeme)
    println("Actual total population size ND*NDeme = ", ND*NDeme)
    println("Area = ",round(Area; sigdigits = 2)," km^2")
    println("Total area from Ademe*ND = ",round(Ademe*ND; sigdigits = 3)," km^2")



    if NDeme==0
        return NaN, NaN, NaN, NaN
    end

    numbuniquetimepoints = length(uniquesamplingtime)


    ##Initialise arrays
    Xt = [fill(NaN + NaN*im, KK) for iy=1:ny, ix=1:nx] #current generation

    

    ncomplexhap = [ [] for ii in eachindex(Nsamp)] #array of an array to hold the complex haplotype data from sampling from simulations, 
                                                    #where Nsamp is an array which holds the number sampled at each distinct location 
                                                    #Note that not all of these will be used as they fall outside the species range


    println("size(ncomplexhap) = ", size(ncomplexhap))

    #initialise demes with WT with frequency 1 in each deme
    for i in eachindex(CartInd) # populate deme cell array with initial wt-frequency of 1 in generation t=1

        # Get linearindex
        k = LinearIndices(CartInd)[i]

        if DemesInAfrica[k] #only initialise demes with Africa shape file (lon,lat) — those outside will always be empty


            Xt[i][1] = 0.0 + 1.0*im

            for kk in eachindex(ncomplexhap)
                ncomplexhap[kk] = Complex[] 
            end


        end


    end

    if !isempty(xmut_target)
        #Copy initialised ncomplexhap
        ncomplexhap_initial = deepcopy(ncomplexhap)
    end

    t = 1 #t is a counter for the generations
    NextGeneration = true #Parameter used to stop the simulation when

    Xsamplemut = 0




    while NextGeneration == true  #Number of generations

        t = t + 1 #Increases the generation counter each iteration

        
        println("t=",t)

        M = [fill(NaN + NaN*im,3*K[i]) for i in CartInd]

     

        truexwt = 0

        if !isempty(xmut_target)
            nmut = 0 #resets every generation
            ncomplexhap = deepcopy(ncomplexhap_initial)
        end

        # @time begin println("1st loop over demes for migration")
        #migration
        for i in CartInd

            # Calculate linear index of CartesianIndex i
            # Get linearindex
            k = LinearIndices(CartInd)[i]
            
            #add if statement that immediately skips this deme if DemesInAfrica (boolean mask) ==0 for this deme
            if DemesInAfrica[k]


                
                hmarker = real(Xt[i][:])
                freq = imag(Xt[i][:])


                # #get x and y indices from CartesianIndex i
                ix = i[2]
                iy = i[1]
                


                indNotNan  = findall(x->!isnan(x),hmarker)


                #loop over all alleles:

                for j in eachindex(indNotNan)
                    #assign each allele to nearest neighbours — 1st element of v below is always the current deme and the rest can be assign to nearest neighbour demes in any order
                    # println("indNotNan[j] =", indNotNan[j])
                    Nj = Int(round(NDeme*freq[indNotNan[j]]))# should change this to enusre that Σj Nj = Ndeme, becuase rounding in this way doesn't!
                    hmarkerj=hmarker[indNotNan[j]]



                    if Nj!=0 #this will effectively remove dead alleles since zero-frequency alleles won't be added to current deme or nearest neighbours


                        if nnn[i]==4
                            v = rand(Multinomial(Nj,[1-p;p/4;p/4;p/4;p/4]));
                        elseif nnn[i]==3
                            v = rand(Multinomial(Nj,[1-3/4*p;p/4;p/4;p/4]));
                        elseif nnn[i]==2
                            v = rand(Multinomial(Nj,[1-1/2*p;p/4;p/4]));
                        elseif nnn[i]==1 #this can rise if nx=1 or ny=1
                            v = rand(Multinomial(Nj,[1-1/4*p;p/4])); 
                        end

                        if sum(v)!=Nj
                            println("v = ",v)
                        end

                        if any(v.<0)
                            println("v = ",v)
                            readline()
                        end


                            
                        #loop over nearest neighbours, separating out wt

                        if v[1]!=0 
                        
                            indFirstNaN = findfirst(x->isnan(x),M[i])

                            if indFirstNaN === nothing #i.e. no NaNs => array full => double size of array for this deme with NaNs

                                M[i] = [M[i];[NaN + im*NaN for hh=1:length(M[i])]]
                                indFirstNaN = findfirst(x->isnan(x),M[i])

                            end

                            M[i][indFirstNaN] = hmarkerj .+ im*v[1] 

                        end

                        
                        

                        
                        for jj in eachindex(nnind[i]) #for the jth allele add migrants calculated in v to nearest neighbour demes

                            if v[jj+1]!=0
                                ii = nnind[i][jj]

                                indFirstNaN = findfirst(x->isnan(x),M[ii])
 

                                if indFirstNaN === nothing

                                    M[ii] = [M[ii];[NaN + im*NaN for hh=1:length(M[ii])]]
                                    indFirstNaN = findfirst(x->isnan(x),M[ii])

                                end


                                M[ii][indFirstNaN] = hmarkerj .+im*v[jj+1]
                            end


                        end
                
                    end

  

                end

                if any(imag.(M[i]).<0)
                    println("M[i] = ",M[i])
                    readline()
                end


        
            end



        end
    
        # end #end for @time 1st migration loop

        # @time begin println("2nd loop over demes for migration")
        #this loop loops though eligible demes to calculate who has migrated into the ith deme
        for i in CartInd #2nd loop over demes for migration


            k = LinearIndices(CartInd)[i]

            if DemesInAfrica[k]

                hmarker = real.(M[i])
                nhmarker = imag.(M[i])



                indNotNan = findall(x->!isnan(x),M[i])

                if !isempty(indNotNan)

                    h = unique(real.(M[i]))
                    h = h[.!isnan.(h)]
                    f = [0.0 for kk=1:length(h)]

                    for j in eachindex(h)

                        if !isnan(h[j])
                            indj = findall(x->x==h[j],hmarker)
                            # sum(nhmarker[indj])
                            f[j] = sum(nhmarker[indj])#/NDeme

                        end


                    end

                    f = f/sum(f[.!isnan.(f)]) #note we are not normalising by NDeme, 
                    #since the population size of the deme due to multinomial migration may have changed. 
                    #We calculate the frequency of allele and assume the deme grows/dies to the carrying 
                    #capacity on a timescale smaller than a generation so the number is NDeme
                    #N.B. if f (& h) only has NaNs this returns NaNs, which means there are nothing in the deme — this is ok given how multinomial sampling works!


                    indNotNanh = findall(x->!isnan(x),h)
                    indNotNanf = findall(x->!isnan(x),f)

                    if length(indNotNanf) != length(indNotNanh)
                        println("h =",h)
                        println("f =",f)
                        readline()
                    end



                    if length(h) > size(Xt[i])[1]
                        # println("size(Xt[i]) = ", size(Xt[i]))
                        # println("Initial array size too small for new migrants: extending array")
                        Xt[i] =  vcat(Xt[i],[NaN + im*NaN for hh=1:K[i]]) #may need to use vcat..
                        # println("size(Xt[i]) = ", size(Xt[i]))

                        K[i]=2*K[i]
                    end



                    if sum(f)>1 +1/NDeme/10
                        println("sum(f) = ",sum(f))
                        println("h = ",h)
                        println("f = ", f)
                        k = LinearIndices(CartInd)[i]             
                        println("Linear index k =",k)
                        println("[iy ix] = ",[i[1],i[2]])
                        

                        println("M[i] =", M[i])

                    end

                    Xt[i][:] = [NaN + im*NaN for hh=1:K[i]] #wipe Xt[i], since the number of alleles can change
                    Xt[i][1:length(h)] = sort(h .+ f*im,by=x->real(x)) 
            
                else #population size after migration is zero: should only happen for very small population sizes per deme

                    Xt[i][:] = [NaN + im*NaN for hh=1:K[i]]
                end


            end


        end

        # end #end of begin for @time

        # @time begin println("loop over demes for drift, selection + mutation")

        numhmarker=0


        for i in CartInd #drift, selection + mutation

            k = LinearIndices(CartInd)[i]

            if DemesInAfrica[k]
                # hmarker = real(Xt[i,:])
                # freq = imag(Xt[i,:])
                hmarker = real(Xt[i][:])
                freq = imag(Xt[i][:])

                indNotNanh = findall(x->!isnan(x),hmarker)
                indNotNanf = findall(x->!isnan(x),freq)

                #if these are empty then this means the population size of deme is 0 and so we can't do anything with this deme until something migrates back in

                # if length(indNotNanf) != length(indNotNanh)
                #     println("hmarker =",hmarker)
                #     println("freq =",freq)
                #     readline()
                # end

                if !isempty(indNotNanh)

              
        
                    indlast = findfirst(x->isnan(x),hmarker)
                    if indlast === nothing #i.e. there are no NaNs in the array - then extend array to allow more alleles
                        
                        freq = [freq;[NaN for hh=1:length(freq)]]
                        hmarker = [hmarker;[NaN for hh=1:length(hmarker)]]
                        indlast = findfirst(x->isnan(x),freq)
                        

                        #also extend size of current demes Xt
                        # println("size(Xt[i]) = ", size(Xt[i]))
                        # println("Initial array size too small for new mutants: extending array")
                        Xt[i] =  vcat(Xt[i],[NaN + im*NaN for hh=1:K[i]])
                        # println("size(Xt[i]) = ", size(Xt[i]))
                        K[i]=2*K[i]


                    end

                    
                    Ki = Int(indlast-1) #last non-NaN value and number of alleles in ith deme

                    indNotNanh = findall(x->!isnan(x),hmarker)
                    indNotNanf = findall(x->!isnan(x),freq)

                    if length(indNotNanf) != length(indNotNanh)
                        println("hmarker =",hmarker)
                        println("freq =",freq)
                        readline()
                    end
                    

                    if isempty(findall(x->x==0.0,hmarker)) #can't check if wt frequency is zero, since it may not be listed in deme any more
                        xwt=0
                        Xi = freq[1:Ki]
                    else
                        xwt = freq[1] 
                        Xi = freq[2:Ki]
                    end
                        
         

                    if 0 < xwt < 1 #then Ki must include the WT
                        # ss = xwt * s
                        ss = xwt * s * (hdom .+ (1 - 2*hdom)*(1-xwt))
                        
                        xx = @. Xi + ss * Xi #Each de novo mutant is given a selection advantage — for xwt=0 ss=0 and the mean frequency doesn't change

                        
                        xx[xx.<0].=0.0
                        xxK = 1 - sum(xx)

                        if xxK < 0 # in case wild type frequency is negative
                            xxK = 0 # set wt frequency to 0
                            xx = xx / sum(xx) #rescales all mutant frequencies to sum to 1
                        end

                        

                        ##Genetic drift
                        if Gaussian == 0 #Sample NDeme number of individuals based on drift, using multinomial sampling
                            z = rand(Multinomial(Int(NDeme), [xxK; xx]))/NDeme
                            
                        else # Sample NDeme number of individuals based on Gaussian approximation of multinomial sampling
                            z = GaussPoissonHybrid_mnrnd(Int(NDeme), [xxK; xx])/NDeme
                            
                        end

                        
                        freq[1:Ki] = z #update frequency vector with values after WF sampling and selection — values after Ki will be NaN
       

                    elseif xwt==0 #deme only has mutants and there can only be genetic drift

                        if !isprobvec(Xi)
                            println("")
                            println("Deme i = ",i)
                            println("Xt[i] = ",Xt[i][:])
                            println("Ki =", Ki)
                            # println("X[i][:,t-1] = ",Xt[i][:])
                            println("M[i] = ",M[i])
                            println("Xi = ",Xi)
                            println("sum(Xi) =",sum(Xi))
                            println("hmarker = ",hmarker)
                            println("freq =",freq)
                            readline()
                        end

                        ##Genetic drift
                        if Gaussian == 0 #Sample NDeme number of individuals based on drift, using multinomial sampling
                            z = rand(Multinomial(Int(NDeme), Xi))/NDeme
                            
                        else # Sample NDeme number of individuals based on Gaussian approximation of multinomial sampling
                            z = GaussPoissonHybrid_mnrnd(Int(NDeme), Xi)/NDeme
                            
                        end

                        freq[1:Ki] = z




                    end 
        
                    

                    #Determining the number of de novo mutants
                    #Note I'm using xwt before WF sampling to determine number of mutants, even though after sampling it could be that xwt=0
                    if xwt > 0 
                        m = pois_rand(NDeme * μ * xwt) #number of de novo mutants generated in each deme (this is based on wild type before selection and before genetic drift)


                        if Ki+m> K[i] # Initial array size too small for new mutants: extending array
                            
                           
                            Xt[i] = vcat(Xt[i],[NaN + im*NaN for hh=1:K[i]])
                            K[i]=2*K[i]

                            freq = [freq;[NaN for hh=1:length(freq)]]
                            hmarker = [hmarker;[NaN for hh=1:length(hmarker)]]

                            # println("i = ",i)
                            # println("K[i] = ",K[i])

                            
                        end


                        if m>0

                            freq[Ki+1:Ki+m] = 1/NDeme*ones(m)
                            hmarker[Ki+1:Ki+m] = collect(mutant_counter +1:mutant_counter+m)
                            
                            #after mutation normalise to frequency
                            freq = freq/sum(freq[.!isnan.(freq)])
                                                         

                            mutant_counter = mutant_counter + m
                        end

                    end

                else

                    xwt=0  #And hmarker and freq will both only have NaNs indicating a dead population in this deme 
                            #— N.B. doing this is not strictly a WF model which would have a fixed population size per deme, 
                            #but this is only really an important distinction if the carrying capacity per deme is very small (NDeme<10)  
                    
                end

                
                Xt[i][:] = hmarker .+ im*freq
               


    
                truexwt = truexwt + xwt

                #List of haplotype markers in this deme
                h0 = hmarker[.!isnan.(freq)]
                #And their corresponding frequencies
                f0 = freq[.!isnan.(freq)] 
                   

                numhmarker = max(numhmarker,length(h0)-1)
                


                if isempty(xmut_target)
                                    # #New Sampling approach

                    # #is current generation equal to one of the times in uniquesamplingtime
                    # indt = findall(x->x.==t,uniquesamplingtime)
                    indt = findall(x->x.==t,samplingtime)

                    if !isempty(indt) #indt is an index to all spacetime lcoations that are sampled in this generation t
                        #now spacetimesample_deme_number(indt) is a list of demes corresponding to sample locations sampled at time t, so find if any equal current deme with linear index k 

                        #create reduced arrays and do all indexing with reference to these
                        Nsampkt  = Nsampt[indt] #these are the number of chromosomes sampled at the space-time sample locations at this time
                        unqspatindx = uniquespatialindx[indt] #these are the corresponding indices to the unique spatial (not space-time) locations
                        indk = findall(x->x.==k,spacetimesample_deme_number[indt]) #note that the indk indices are into the reduced arrays and indk may have more than one sample site since there may be mulitple sites per deme

                        #Now loop over these 

                        for kk in indk
                            Ksampt = Int(only(Nsampkt[kk]))
                            nnt = rand(Multinomial(Ksampt,f0)) #Sample the haplotype frequencies with Ksampt, if Ksampt =0, then nnt=0
                                
                            for hi in eachindex(h0)
                                push!(ncomplexhap[unqspatindx[kk]], h0[hi] .+ im*nnt[hi]) #this adds contributions to the same sample point in the ith deme

                            end

                        end


                    end
                else #xmut_target is not empty and change sampling method
                    indk = findall(x->x.==k,spacetimesample_deme_number)
                    

                    if !isempty(indk) # i.e. there is at least one space-time sample site which matches the current deme with linear index k


                        #Now loop over these space-time sites

                        for kk in indk
                            # println("kk = ",kk)
                            Ksampt = Int(only(Nsampt[kk]))

                            nnt = rand(Multinomial(Ksampt,f0)) #Sample the haplotype frequencies with Ksampt, if Ksampt =0, then nnt=0
                                
                            for hi in eachindex(h0)
                                push!(ncomplexhap[uniquespatialindx[kk]], h0[hi] .+ im*nnt[hi]) #this adds contributions to the same sample point in the ith deme

                                if h0[hi]!=0.0 # => mutant haplotype
                                    # println("hello")
                                    nmut = nmut + nnt[hi] #this accumulates all mutant haplotypes sampled over sample points and times
                                end
                                # println("nmut = ",nmut)
                            end

                        end
                    end
                end
                    

                

            end
        end

        # end #end for begin block for @time

        truexwt = truexwt/ND

        if !isempty(xmut_target)
            Xsamplemut = nmut/sum(Nsampt)
        end
        

        truefreq = push!(truefreq,1-truexwt)

        

        if !isempty(xmut_target)
            #Temporary for production of figures for paper — comment out normally and uncomment t>=maxt condition
            if Xsamplemut >= xmut_target
                NextGeneration = false
            end
        else

            
            if t>=maxt
                NextGeneration = false
            end

        end

    end #while NextGeneration == true

    println("t=",t)
    # println("Sample WT freq=", XWT)
    println("Sample mutant freq=", Xsamplemut)

   



    return Xt, ncomplexhap, nx, Ademe, NDeme, t, truefreq


end #function





function CalculateDemes(Ademe, Area, Ratio, lon_shape,lat_shape,R)



    if !isempty(lon_shape) #Area not used

    
        x_shape,y_shape = LongLat2km(lon_shape,lat_shape,R)

        minx = minimum(x_shape)
        miny = minimum(y_shape)

        # #Simulations performed with offset in x and y 
        x_shape = x_shape .- minx
        y_shape = y_shape .- miny

        maxx = maximum(x_shape)
        maxy = maximum(y_shape)

        Area = maxx*maxy; #not actual area of simulation but smallest rectangle that would fit
        Ar = maxx/maxy;

        println("Longitudinal and latitudinal shape file input. Based on these:")
        println("Δx = ", round(maxx;digits=1)," km")
        println("Δy = ", round(maxy; digits=1)," km")
        println("Area calculated from these Δx*Δy = Area =", round(Area;sigdigits=3)," km^2")
        println("Aspect ratio Δx/Δy = ", round(Ar;sigdigits=3))

        NNDD = Area/Ademe

        #nnx and nny are not the final number of demes in x and y direction, but they correspond to the number of demes *if* the simulation area just encapsulates the shape file
        #however, for convenience I want to choose deme centres to be anywhere within or on the edge of the shape file, and so there will be a proportion of area outside the shape file boundary
        #but nx = nnx+1 and ny = nny+1
        nny = round(Int,sqrt(NNDD/Ar))
        nnx = round(Int,Ar*nny)

        println("nnx=",nnx)
        println("nny=",nny)
        # NDD = nnx*nny
        println("Total number of demes from shape file after correcting for Aspect Ratio = ",nnx*nny)
        AA = Ademe*nnx*nny

        
        #Now, demes centered on exactly or near edge of shape outline in lon, lat, will have a proportion of their area strictly outside the shape file & => we need to calculate the maximum extent of that area A and the number of  
        A = AA + Ademe*(nnx+nny+1) #actual area to be simulated:

        nx = nnx + 1
        ny = nny + 1

        println("nx=",nx)
        println("ny=",ny)

        ND = nx*ny
        println("Total number of demes in simulation ND = ",ND)
        
        println("Modified total area for simulations (km^2) =",round(A; sigdigits = 3))
        println("Rectangular area just encompassing shape file (km^2) = ", round(Area; sigdigits = 3))

        println("Modified aspect ratio Ar = nx/ny = ",round(nx/ny;sigdigits = 3))
        println("Aspect ratio of rectangular area just encompassing shape file = ",round(Ar;sigdigits=3))

    else

        println("Aspect Ratio Δx/Δy = ", Ratio)

        NND = Area/Ademe; 
        println("1st estimate of number of demes: Area/Ademe =",NND)
        ny = round(Int,sqrt(NND/Ratio))
        nx = round(Int,Ratio*ny)

        println("nx=",nx)
        println("ny=",ny)
        ND = nx*ny
        println("Total number of demes after correcting for specified Aspect Ratio = ",ND)
        A = Ademe*ND
        println("Modified total area (km^2) =",round(A; sigdigits = 3))
        println("Requested area (km^2) = ", round(Area; sigdigits = 3))

        println("Modified aspect ratio Ar = nx/ny = ",round(nx/ny;sigdigits=3))
        println("Requested aspect ratio Ar = ",Ratio)

        x_shape = NaN
        y_shape = NaN 
        minx = NaN
        miny = NaN



    end



    return ND, nx, ny, x_shape, y_shape, minx, miny
end





function CalcNearestNeighbours(nx,ny,DemesInAfrica)

    #Calculate array of an array of nearest neighbour indices for demes & update DemesInAfrica for any demes with 0 nearest neighbours
    #e.g. if they are on a peninsula
    nnind = [fill(Int(0),4) for iy=1:ny, ix=1:nx] #current generation
    #arrays of number of nearest neighbours for each deme
    nnn = zeros(ny,nx)

    CartInd = CartesianIndices((1:ny,1:nx))

    for i in CartInd

        # #get x and y indices from CartesianIndex i
        ix = i[2]
        iy = i[1]

        # println("Subscript indices [iy,ix] = ",[iy,ix])

        # Get linearindex
        k = LinearIndices(CartInd)[i]
  

        #nearest neighbour index
        nnsubind=[[iy,ix-1],[iy-1,ix],[iy+1,ix],[iy,ix+1]]
        # println("nnsubind = ",nnsubind)
        nnindi = fill(NaN,1,4) #nearest neighbour indices
        #check bounds and convert to linear/CartesianIndex
        for kk=1:4
            iix = nnsubind[kk][2]
            iiy = nnsubind[kk][1]

            iik = iiy+(iix-1)*ny

            # println("iix = ",iix)
            # println("iiy = ",iiy)

            

            if (0 < iiy <=ny) & (0 < iix <= nx) #*** once mapping onto sub-Saharan Africa will need to modify this to include excluded demes
                #e.g. with a boolean mask/array called say "in" which is 1 if within shape file and 0 if outside 
                # println("hello")


                #if this deme is within the rectangular grid/bounding box then need to check it is also within Africa
                #can't do it in the reverse order since if a migratory deme is outside the rectangular 
                #grid then the linear index will correspond to nonsense wrt to the linear indexing in DemesInAfrica 
                if DemesInAfrica[iik]

                    nnindi[kk] = iiy + (iix-1)*ny  

                end 
                        
            end
       
            

        end

        # nnind = nnind[findall(x->!isnan(x),nnind)] #simple indexing nnind=nnind[.!isnan(nnind)] doesn't work!!
        nnind[i] = Int.(nnindi[.!isnan.(nnindi)])
        # nnind[i] = Int.(nnind[i]) #this should now contain nearest neighbour indices inc current in linear index form


        if !isempty(nnind[i])
            # println("Nearest neighbour indices = ",nnind)
            nnn[i] = length(nnind[i])
        else #if empty the deme whilst in Africa shape, has no nearest neighbour demes in shape (i.e. on coast and oddly jutting out peninsula)
            DemesInAfrica[k] = false

        end

        



    end

    return nnind, nnn, DemesInAfrica


end



function LongLat2km(lon,lat,R)

    #lon is the longitudinal angle from meridian in degrees
    #latitude is the angle up or down from equator -90<=lat<=90 (degrees), where -90 is the south pole and 90 the north
    #θ is the spherical polar coordinate angle down from vertical 0<=θ<=π, where θ=0 is the north pole, and θ=π is the south pole
    #need to convert lat->θ: => lat' = 90-lat — **Don't need to do this**
    φ = pi/180*lon
    θ = pi/180*lat
    # θ = pi/180*(90 .- lat)
    
    
    # % %x-distance (assume latitude fixed and calculate distance to 0-longitude)
    # % psi = acos(cos(φ).*cos(θ).^2 + sin(θ).^2);
    # % psi = acos((cos(φ)-1).*cos(θ).^2+1);
    
    
    # %Havesine version: meant to be more numerically accurate:
    # %hav(ψ) = hav(Δθ) + cos(θ1)cos(θ2)*hav(Δφ)
    # %fix θ=>Δθ=0 and hav(0) = 0 => hav(ψ) = cos(θ)^2*hav(Δφ)
    
    ψ = invhav.(cos.(θ).^2 .*hav.(φ))
    
    x = R*ψ.*sign.(φ)
    
    # %y-distance (assume fixed longitude: psi=θ)
    y = R*θ
    
    return x, y
end
    
    
function hav(θ)

    h = sin.(θ/2).^2;

return h

end
    
function invhav(x)

    # %sqrt(h) = sin(θ/2) 
    # %=> 2*asin(sqrt(h) = θ

    θ = 2*asin.(sqrt.(x));

    return θ

end