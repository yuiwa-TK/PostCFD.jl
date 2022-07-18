using DelimitedFiles
"""
read a file of "header&data format" written in "little-endian" & "stream".  
this function returns Tuple(labels::strings, digit:Num)
"""
function read_header(filename::String)
    # HEADER があるfileを読み，数値データを配列で返す．
    data = readdlm(filename)
    labels = data[1,:]
    digits = data[2:end,:]
    return labels,digits
end

