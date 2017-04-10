function [dist] = fhindeptest3(x,y,num_shuffle)

n = length(x);

UV = [pobs(x) pobs(y)];

dist = zeros(n,(num_shuffle+1)*4);

dist(:,1) = point_to_line(UV, repmat([0,0],n,1), repmat([1,1],n,1));
dist(:,2) = point_to_line(UV, repmat([0,1],n,1), repmat([1,0],n,1));   % could also be done w/ a fliplr
dist(:,3) = point_to_line(UV, repmat([0.1,0],n,1), repmat([0.9,1],n,1));
dist(:,4) = point_to_line(UV, repmat([0,0.9],n,1), repmat([1,0.1],n,1));


for ii=1:num_shuffle
    UV_shuffled = shuffle_pobs(UV);
    dist(:,ii*4+1) = point_to_line(UV_shuffled, repmat([0,0],n,1), repmat([1,1],n,1));
    dist(:,ii*4+2) = point_to_line(UV_shuffled, repmat([0,1],n,1), repmat([1,0],n,1));
    dist(:,ii*4+3) = point_to_line(UV_shuffled, repmat([0.1,0],n,1), repmat([0.9,1],n,1));
    dist(:,ii*4+4) = point_to_line(UV_shuffled, repmat([0,0.9],n,1), repmat([1,0.1],n,1));
end

end

