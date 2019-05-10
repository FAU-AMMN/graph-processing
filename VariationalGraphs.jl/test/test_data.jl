using DelimitedFiles
#------------------------------------------------
tmp = readdlm("data.txt")
data = Array{Float64,3}(undef, 4, 4, 3)
for i = 1:3
    data[:, :, i] = tmp[((i - 1) * 4 + 1):(i * 4), :]
end
    
