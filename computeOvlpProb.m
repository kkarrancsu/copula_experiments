function [ovlp] = computeOvlpProb(mu1,sigma1,mu2,sigma2)

polyroots = getNormPDFRoots(mu1,mu2,sigma1,sigma2);
if(isempty(polyroots))
    ovlp = 0;
elif(len(polyroots)==1)
    ovlp = 1-normcdf(polyroots(1),mu1,sigma1)+normcdf(polyroots(1),mu2,sigma2);
else
    if(polyroots(1)>polyroots(2))
        c2 = polyroots(1);
        c1 = polyroots(2);
    else
        c2 = polyroots(2);
        c1 = polyroots(1);
    end
    if(sigma2>sigma1)
        F2mu = mu2; F2sigma = sigma2;  
        F1mu = mu1; F1sigma = sigma1;
    else
        F2mu = mu1; F2sigma = sigma1;
        F1mu = mu2; F1sigma = sigma2;
    end
    ovlp = normcdf(c2,F2mu,F2sigma) - ...
           normcdf(c1,F2mu,F2sigma) + ...
           normcdf(c1,F1mu,F1sigma) + ...
           (1- normcdf(c2,F1mu,F1sigma));
end

end

function [polyroots] = getNormPDFRoots(m1,m2,std1,std2)
    a = 1/(2*std1^2) - 1/(2*std2^2)
    b = m2/(std2^2) - m1/(std1^2)
    c = m1^2 /(2*std1^2) - m2^2 / (2*std2^2) - log(std2/std1)
    polyroots = roots([a,b,c])
end