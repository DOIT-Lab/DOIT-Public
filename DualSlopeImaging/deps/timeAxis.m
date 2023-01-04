function [ t ] = timeAxis( fs, l )
%[ t ] = timeAxis( fs, l )

    t=(0:(1/fs):((l-1)/fs))';
end

