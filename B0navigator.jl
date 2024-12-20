export eddyCurrentCorrection, deltaB0EstimationLL, deltaB0Estimation_iterativeFitting

using PyCall
using LsqFit
using FFTW
using LinearAlgebra
using StatsBase
using NaNStatistics

function eddyCurrentCorrection(data::Array{Complex{T}}, uniqueKy::Vector{Int}) where T<:AbstractFloat
    @debug("Correct AC data for hardware imperfections.")
    ## Expected size of AC region as input: size(data) = (numSpokes, numKz, numEchoes, numTI, numCoils)
    ## uniqueKy: Order of spokes based on the angle increment
    shape = size(data)
    numSpokes = shape[1]
    numSlices = shape[2]
    numEchoes = shape[3]
    numTI = shape[4]
    numCoils = shape[5]

    ## IFFT along ky
    data = fftshift(ifft!(ifftshift(data,2),2), 2)

    ## Change to profile order
    ky = uniqueKy
    data[ky.+1,:,:,:,:] = data
    data = permutedims(data, [1,2,3,5,4])
    data = reshape(data, numSpokes, numCoils*numTI*numEchoes*numSlices)

    ## Do correction according to Rosenzweig et al.
    nH = 5
    phi = diff(collect(2*pi.-2*pi.*(0:numSpokes-1)/numSpokes))[1]
    nCorr = ones(ComplexF32, numSpokes, 2*nH)
    for i in 1:2*nH
        for j in 1:numSpokes
            nCorr[j, i] = exp(1im*(-1)^i*ceil(i/2)*(j-1)*phi)
        end
    end
    data = reshape(transpose(data - nCorr*(pinv(nCorr)*data)), numCoils*numEchoes*numSlices, 
        numTI, numSpokes)
    data = permutedims(data, [1, 3, 2])

    ## Change back to temporal order
    data = data[:,ky.+1,:]
    return data
end

function deltaB0EstimationLL(data::Array{Complex{T}}, params::Dict{Symbol,Any}; verbose::Bool=true) where T<:AbstractFloat
    ## Pre-processing
    if :uniqueKy in keys(params)
        shape = size(data)
        numSpokes = shape[1]
        numSlices = shape[2]
        numEchoes = shape[3]
        numTI = shape[4]
        numCoils = shape[5]
        data = eddyCurrentCorrection(data, params[:uniqueKy])
    else
        ## Expected size of AC region: size(data) = (numEchoes, numSpokes, numTI) -> numCoils = 1, numSlices = 1
        shape = size(data)
        numSpokes = shape[2]
        numTI = shape[3]
        numEchoes = shape[1]
        numCoils = 1
        numSlices = 1
    end

    @debug("Estimate B0 variations.")
    data = reshape(permutedims(data, [3,2,1]), numSpokes*numTI, :, 1, numEchoes, numCoils) # Reshape for LL acquisition
    numSlices = size(data,2)

    ## Initilize graph-cut parameters
    hmrGC = pyimport("hmrGC.dixon_imaging")
    th = 1.0 
    paramsGandalf = Dict{String, Any}()
    if :fatModel in keys(params)
        paramsGandalf["FatModel"] = params[:fatModel]
    end
    paramsGandalf["TE_s"] = params[:TE_s] # echo times in s
    paramsGandalf["centerFreq_Hz"] = params[:centerFreq_Hz] # center frequency in Hz
    paramsGandalf["fieldStrength_T"] = params[:fieldStrength_T]  # field strength in T
    paramsGandalf["voxelSize_mm"] = [1.5, 1.5, 5.5]  # hard-coded, does not have an effect on B0 navigator

    phasormap = zeros(Float64, size(data, 1), size(data, 2), size(data, 3), size(data, 5))
    water = zeros(ComplexF64, size(data, 1), size(data, 2), size(data, 3), size(data, 5))
    fat = zeros(ComplexF64, size(data, 1), size(data, 2), size(data, 3), size(data, 5))
    mip = maximum(abs.(data), dims=4)[:,:,:,1,:]
    mask = zeros(Bool, size(mip,1),size(mip,2),size(mip,3),size(mip,4))

    for i in eachindex(axes(mip,4))
        mask[:,:,:,i] = repeat(maximum(mip[:,:,:,i], dims=1) .> (th / 100 * percentile(reshape(maximum(mip, dims=1),:),95)),  
                outer=[size(mip,1), 1, 1, 1]) ## better choose maximum for LL

        if sum(mask[:,:,:,i]) > 0
            ## Graph-cut to solve for B0 estimate
            obj = hmrGC.MultiEcho(data[:,:,:,:,i], mask[:,:,:,i], paramsGandalf)
            obj.range_fm_ppm = params[:range_fm_ppm]
            obj.sampling_stepsize_fm = params[:sampling_stepsize_fm] 
            obj.verbose = false
            obj.perform("single-res")
            phasormap[:,:,:,i] = obj.fieldmap
            water[:,:,:,i] = obj.images["water"]
            if "fat" in keys(obj.images)
                fat[:,:,:,i] = obj.images["fat"]
            end
        end
    end

    output = Dict{Symbol, Any}()
    output[:deltaB0_raw_water] = water
    output[:deltaB0_raw_fat] = fat
    output[:deltaB0_raw_b0] = phasormap
    return output
end

function deltaB0Estimation_iterativeFitting(output::Dict{Symbol,Any}, numContr::Int) 
    @debug("Estimate motion curve from B0 navigator.")
    iterations = 40

    phasormap = deepcopy(output[:deltaB0_raw_b0])
    water = deepcopy(output[:deltaB0_raw_water])
    mask_signal = repeat(sum(abs.(water), dims=1) .== 0, outer=[size(phasormap,1),1,1,1]) ## Mask slices/coils without signal
    phasormap[mask_signal] .= NaN
    indices = zeros(Int, 2, size(phasormap, 2), size(phasormap, 4))
    indices[1,:,:] = repeat(reshape(1:size(phasormap, 2), :, 1), outer=[1, size(phasormap, 4)])
    indices[2,:,:] = repeat(reshape(1:size(phasormap, 4), 1, :), outer=[size(phasormap,2), 1])
    indices = reshape(indices, 2, :) # Indices to identify slice and coil number 

    numSlices = size(phasormap,2)
    numSpokes = Int(size(phasormap,1) ./ numContr)
    samples = reshape(phasormap, size(phasormap,1), numSlices*size(phasormap, 4))

    samples -= repeat(samples[end:end, :], outer=[size(samples,1),1]) # Remove fm offset for last contrast and acquired shot
    motionSignal = zeros(size(samples,1))  # Initilize motionSignal
    t = collect(range(0, stop=1, length=size(samples,1)))
    t = t.-1

    params = nothing
    for iteration in 1:iterations
        @debug("Iteration: ", iteration)
        if params == nothing
            params = zeros(2+numContr, size(samples,2))
        end
        resid = zeros(size(samples,2))
        model_sample(t, p) = vec(repeat(p[2:2+numContr-1], outer=[numSpokes,1])) .+ p[1]*t # first model to fit contrast offset and drift 
        model_motion(t, p) = p[1] * motionSignal # second model to fit B0 variation amplitude
        for i in 1:size(samples,2)
            if !any(isnan.(samples[:,i]))
                if iteration == 1
                    fit = curve_fit(model_sample, t, samples[:,i], params[2:end,i]) # first model is only fitted once, independent of motionSignal
                    params[2:end,i] = fit.param
                end
                fit = curve_fit(model_motion, t, samples[:,i] - model_sample(t, params[2:end,i]), params[1:1,i]) # fit B0 variation amplitude
                params[1:1,i] = fit.param[1:1]
                ## Compute residual
                if abs(params[1,i]) > 1
                    resid[i] = 1/abs(params[1,i]) * LsqFit.mse(fit) 
                else
                    resid[i] = LsqFit.mse(fit) 
                end
            end
        end

        mean_resid = mean(resid)
        q3 = nanpctile(resid[:], 75)
        iqr = nanpctile(resid[:], 75) - nanpctile(resid[:], 25)
        std_resid = std(resid)
        # @show mean_resid

        outlier = resid .> q3+1.5*iqr # Threshold to remove outlier 
        # @show outlier
        if any(outlier) || iteration == 1
            ## Remove outliers
            if iteration > 1
                @debug("Outliers found: ", sum(outlier))
                samples = samples[:, .!outlier]
                params = params[:, .!outlier]
                indices = indices[:, .!outlier]
            end

            sample_demodulated = deepcopy(samples)
            for i in 1:size(samples,2) # demodulate fitted contrast variations and B0 drift from samples
                sample_demodulated[:,i] = samples[:,i] .- model_sample(t, params[2:end,i])
            end
            sign_alpha = sign.(params[1:1,:])
            sign_alpha[sign_alpha .== 0] .= 1

            motionSignal = median(sample_demodulated .* sign_alpha, dims=2)[:,1] # compute new estimate for motionSignal
        else
            @debug("No outliers found")
            break
        end
    end

    ## Normalize motion signal 
    amplitude_motionSignal = (quantile(motionSignal, 0.95) - quantile(motionSignal, 0.05)) ./ 2
    motionSignal = motionSignal ./ amplitude_motionSignal    
    params[1:1,:] = params[1:1,:] .* amplitude_motionSignal

    ## Reorder estimated parameters
    params_reordered = zeros(size(params,1), size(phasormap,2), size(phasormap,4))*NaN
    for ind = 1:size(indices,2)
        slice = indices[1,ind]
        coil = indices[2,ind]
        params_reordered[:,slice,coil] = params[:,ind]
    end
    output[:deltaB0_iterativeFitting_motionCurve] = motionSignal
    output[:deltaB0_iterativeFitting_params] = params_reordered
    return output
end