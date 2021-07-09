function strobes = fix_missing_strobes(strobes)
%****** function strobes = fix_missing_probes(strobes)
%****** 
%****** there was an error in our first data file where
%****** codes of 64 were written as 0, so we had to replace
%****** those 0's with 64's, but had to make sure they were
%****** part of the expected start/end probe sequence

    zz = find( strobes == 0);
    zN = size(zz,1);
    N = size(strobes,1);
    for k = 1:zN
        kk = zz(k);
        kN = kk + 7;
        if (kN <= N)   
           if (strobes(kN) == 63)
              strobes(kk) = 64;
           end
        end        
    end

return;