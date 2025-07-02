# Written by Rebekka Muller

using LinearAlgebra
using DataFrames
using CSV

### ---------------- AMINOACIDS, CODONS, and FITNESS ---------------- ###
# amino acids
const aa = [
    "Stop",
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Gln",
    "Glu",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val",
]

# map Φ, dictionary with codons as keys and amino acids as corresponding values
const Φ =
    Dict(1 => "Lys", 2 => "Asn", 3 => "Lys", 4 => "Asn",
        5 => "Thr", 6 => "Thr", 7 => "Thr", 8 => "Thr",
        9 => "Arg", 10 => "Ser", 11 => "Arg", 12 => "Ser",
        13 => "Ile", 14 => "Ile", 15 => "Met", 16 => "Ile",
        17 => "Gln", 18 => "His", 19 => "Gln", 20 => "His",
        21 => "Pro", 22 => "Pro", 23 => "Pro", 24 => "Pro",
        25 => "Arg", 26 => "Arg", 27 => "Arg", 28 => "Arg",
        29 => "Leu", 30 => "Leu", 31 => "Leu", 32 => "Leu",
        33 => "Glu", 34 => "Asp", 35 => "Glu", 36 => "Asp",
        37 => "Ala", 38 => "Ala", 39 => "Ala", 40 => "Ala",
        41 => "Gly", 42 => "Gly", 43 => "Gly", 44 => "Gly",
        45 => "Val", 46 => "Val", 47 => "Val", 48 => "Val",
        49 => "Stop", 50 => "Tyr", 51 => "Stop", 52 => "Tyr",
        53 => "Ser", 54 => "Ser", 55 => "Ser", 56 => "Ser",
        57 => "Stop", 58 => "Cys", 59 => "Trp", 60 => "Cys",
        61 => "Leu", 62 => "Phe", 63 => "Leu", 64 => "Phe"
        )

# dictionary with amino acids as keys and codons as corresponding values - the inverse of Φ
const aminoacids = Dict(a => sort!(findall(==(a), Φ)) for a in unique(values(aa)))

# definition of the fixation rate ω_γ
function ω(γ)
    if γ == 0
        ω = 1.0
    else
        ω = γ / (1 - exp(-γ))
    end
end

### ---------------- MARKOV PROCESS ---------------- ###

"""
    SingleSite(α, μ, aminoacid_fitness)

Computes the stationary distribution of a codon frequencies at a single site for a given amino acid fitness landscape.

## Arguments

- α:                    parameter for the discrete MC
- μ:                    mutation intensity
- aminoacid_fitness:    scaled Malthusian amino acid fitness landscape for one codon site
"""
function SingleSite(α, μ, aminoacid_fitness)

    if length(aminoacid_fitness) != 21
        @info("Fitness vector for amino acids has the wrong length.")
    end

    ### ----------- CONSTRUCTION OF TRANSITION MATRICES ----------- ###

    β = (1 - α) / 2

    # construct H (transition probability matrix for nucleotide changes)
    H = [0 β α β; # A
         β 0 β α; # C
         α β 0 β; # G
         β α β 0] # T

    # construct S 
    Id4 = Matrix{Float64}(I, 4, 4)
    S = [H β*Id4 α*Id4 β*Id4;
         β*Id4 H β*Id4 α*Id4;
         α*Id4 β*Id4 H β*Id4;
         β*Id4 α*Id4 β*Id4 H]

    # construct M (transition probability matrix for codon changes)
    Id16 = Matrix{Float64}(I, 16, 16)
    M = 1/3 .* [S β*Id16 α*Id16 β*Id16;
                β*Id16 S β*Id16 α*Id16;
                α*Id16 β*Id16 S β*Id16;
                β*Id16 α*Id16 β*Id16 S]

    # check doubly stochasticity of M
    #@test all(x -> sum(x) ≈ 1, eachrow(M)) # should be ones everywhere: println("Rows of M sum to one.")
    #@test all(x -> sum(x) ≈ 1, eachcol(M)) # should be ones everywhere: println("Columns of M sum to one. Hence M is doubly stochastic.")

    # decomposition
    # synonymous mutations between sense codons
    M_syn = [
        if u != v && Φ[u] == Φ[v] && Φ[v] != "Stop"
            M[u, v]
        else
            0.0
        end
        for u in 1:64, v in 1:64
        ]

    # nonsynonymous mutations between sense codons
    M_non = [
        if u != v && Φ[u] != Φ[v] && Φ[v] != "Stop" && Φ[u] != "Stop"
            M[u, v]
        else
            0.0
        end
        for u in 1:64, v in 1:64
        ]
    
    # mutations involving stop codons
    M_stop = [
        if u != v && (Φ[v] == "Stop" || Φ[u] == "Stop")
            M[u, v]
        else
            0.0
        end
        for u in 1:64, v in 1:64
    ]

    #@show findall(!iszero,M_stop)
#=
    @test M == M_syn + M_non + M_stop   # println("It holds M = M_syn + M_non + M_stop.")
    @test issymmetric(M)
    @test issymmetric(M_syn)       
    # println("M and M_syn are symmetric. M_non and M_stop are not symmetric by construction.")     

    @test iszero(diag(M))
    @test iszero(diag(M_syn))
    @test iszero(diag(M_non))
    @test iszero(diag(M_stop))
=#

    ### ----------- CONSTRUCTION Ω ----------- ###

    Ω = Matrix{Float64}(undef, 64, 64)
    for (i, u) in enumerate(aa)
        u_codons = aminoacids[u]
        if u !== "Stop"
            # Mutation from a sense codon
            F_u = aminoacid_fitness[i]
            for (j, v) in enumerate(aa)
                # Compute fixation rate (scaled by N)
                ω_γ = if v == u
                    # Synonymous mutation
                    1.0
                elseif v != "Stop"
                    # Non-synonymous mutation to a sense codon
                    F_v = aminoacid_fitness[j]
                    ω(F_v - F_u)
                else
                    # Non-synonymous mutation to a stop codon
                    0.0
                end

                v_codons = aminoacids[v]
                for b in v_codons, a in u_codons
                    Ω[a, b] = ω_γ
                end
            end
        else
            # Mutation from a stop codon
            for v in aa
                # Compute fixation rate (scaled by N)
                ω_γ = if v == u
                    # Synonymous mutation
                    1.0
                else
                    # Non-synonymous mutation to a sense codon
                    0.0
                end

                v_codons = aminoacids[v]
                for b in v_codons, a in u_codons
                    Ω[a, b] = ω_γ
                end
            end
        end
    end

    #@test all(!isnan, Ω)

    ### ----------- CONSTRUCTION OF MARKOV PROCESS ----------- ###

    # infinitesimal generator
    V_diag = Diagonal(vec(sum(Ω .* M; dims = 2)))
    Q = 3 .* μ .* (Ω .* M .- V_diag)
    #@test all(iszero, diag(Ω .* M))

    # Is rowsum zero for Q?
    #@test all(x -> isapprox(sum(x), 0.0; atol = 10 * eps()), eachrow(Q))    #println("Rows of Q sum to zero.")

    ### ----------- STATIONARY DISTRIBUTION ----------- ###

    sense_codons = [i in aminoacids["Stop"] ? 0.0 : 1.0 for i in 1:64]
    stop_codons = [i in aminoacids["Stop"] ? 1.0 : 0.0 for i in 1:64]
    η_sense = vcat(Q', sense_codons', stop_codons') \ vcat(zeros(64), 1.0, 0.0)     # set weights on stop codons to zero
    η_stop = vcat(Q', sense_codons', stop_codons') \ vcat(zeros(64), 0.0, 1.0)      # set weights on sense codons to zero

    #=
    @test isapprox(sum(η_sense), 1.0, atol = 10 * eps())  #println("The constructed vector η is a proability distribution.")
    @test maximum(abs, Q' * η_sense) ≈ 0.0 atol=1e-9 
    @test isapprox(sum(η_stop), 1.0, atol = 10 * eps())  #println("The constructed vector η is a proability distribution.")
    @test maximum(abs, Q' * η_stop) ≈ 0.0 atol=1e-9 

    @test isapprox(η_sense[aminoacids["Stop"]], zeros(3), atol = 10*eps())          #println("Stationary weights on stop codons are zero.")
    #η_sense[aminoacids["Stop"]] .= 0.0
    =#

    ### ----------- dN/dS ----------- ###

    ηt = transpose(η_sense)
    Q_non = @. 3 * μ * (Ω * M_non)                  # rates for nonsynonymous substitutions
    Q_syn = @. 3 * μ * (Ω * M_syn)                  # rates for synonymous sustitutions
    KN = sum(ηt * Q_non)                            # now D_N in the manuscript
    KS = sum(ηt * Q_syn)                            # now D_S in the manuscript
    LN = sum(ηt * (3 .* μ .* M_non))                # now δ_N in the manuscript
    LS = sum(ηt * (3 .* μ .* M_syn))                # now δ_S in the manuscript
    dNdS = (LS / LN) * KN / KS

    #@test(isapprox(KS,LS))

    return (; η_sense, η_stop, dNdS, KN, KS, LN, LS)
end

#= η_sense = SingleSite(0.5, 1e-6, aafitness(1e-3, 1e-3, "Ala"))[1]
η_sense[aminoacids["Stop"]]
η_sense[aminoacids["Ala"]]

η_stop = SingleSite(0.5, 1e-6, aafitness(1e-3, 1e-3, "Ala"))[2]
η_stop[aminoacids["Stop"]]
η_stop[aminoacids["Ala"]]  =#

const CODONS = [z * y * x for z in "ACGT" for y in "ACGT" for x in "ACGT"]

struct Wrightian
    Nₕ::Int
end

struct Malthusian end

const FitnessScaling = Union{Wrightian,Malthusian}

"""
    stationary_distribution(α, μ, FitnessLandscape, scaling::FitnessScaling)

Computes the stationary distribution of a codon sequence for a given amino acid fitness landscape.

## Arguments

- α:                    parameter for the discrete MC
- μ:                    mutation intensity
- FitnessLandscape:     amino acid fitness landscape for one codon sequence
- scaling:              if `FitnessLandscape` is scaled Malthusian `Malthusian()`, if `FitnessLandscape` is Wrightian `Wrightian(Nₕ)`
"""
function stationary_distribution(α, μ, FitnessLandscape, scaling::FitnessScaling)
    # read in fitness landscape
    fitness = CSV.read(FitnessLandscape, DataFrame; delim = ',')

    if scaling isa Wrightian
        # convert Wrightian fitness to scaled Malthusian fitness
        transform!(fitness, Not(["position", "prefAA"]) .=> ByRow(x -> 2 * scaling.Nₕ * log(x)); renamecols = false)
    end

    # output data frame, same columns for position and prefAA as in input df
    # add column for each codon with stationary weights and columns for dNdS-quantities
    sequence = select(fitness, :position, :prefAA, AsTable(aa) => ByRow(x -> (res = SingleSite(α, μ, x); vcat(res.η_sense, res.dNdS, res.KN, res.KS, res.LN, res.LS))) => vcat(CODONS, "dNdS", "KN", "KS", "LN", "LS"); copycols = false) 

    stdistr = select(sequence, Not("dNdS", "KN", "KS", "LN", "LS"))
    dNdS = select(sequence, :position, :prefAA, "dNdS", "KN", "KS", "LN", "LS")
    FileName = replace("$(FitnessLandscape)", ".fl" => ".sd")
    # Save table to disk
    CSV.write(joinpath(@__DIR__, "$(FileName)"), stdistr)
    CSV.write(joinpath(@__DIR__, "dNdS_for_$(FitnessLandscape)"), dNdS)

    # get sequence dNdS
    sequence_dNdS = (sum(dNdS[!, "LS"]) / sum(dNdS[!, "LN"])) * (sum(dNdS[!, "KN"]) / sum(dNdS[!, "KS"]))
    
    return stdistr, dNdS, sequence_dNdS
end

#stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.001_L=1000_Nh=500.csv", Malthusian())
#stationary_distribution(0.5, 1e-6, "WrightianFitnessLandscape_s=0.001_L=1000.csv", Wrightian(500))
stationary_distribution(0.5, 1e-6, "AA_Fitness_Landscape.fl", Wrightian(5e2))



#=
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.0001_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.0005_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.001_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.0025_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_s=0.00125_L=1000_Nh=500.csv")


stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_Uniform_smax=0.0001_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_Uniform_smax=0.0005_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_Uniform_smax=0.001_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_Uniform_smax=0.0025_L=1000_Nh=500.csv")
stationary_distribution(0.5, 1e-6, "ScaledMalthusianFitnessLandscape_Uniform_smax=0.005_L=1000_Nh=500.csv")

stationary_distribution(0.5, 1e-6, "ScaledMalthusianInfluenzaFitness.csv")

=#

