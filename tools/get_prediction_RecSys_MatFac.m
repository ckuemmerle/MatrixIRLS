function predict = get_prediction_RecSys_MatFac(U,V,S,rowind_test, ...
    colind_test,centscale,clip_ratingscale)

m=length(rowind_test);
data=struct;
data.rowind = rowind_test;
data.colind = colind_test;
data.d1     = size(U,1);
data.d2     = size(V,1);

U = U*S;
predict =  partXY(U', V', rowind_test, colind_test, m)';

if not(isempty(centscale))
    predict = centscale_data(predict,data,centscale,'backward');
end

if not(isempty(clip_ratingscale))
    if length(clip_ratingscale) > 1
        ratingscale_min = clip_ratingscale(1);
        ratingscale_max = clip_ratingscale(2);
        predict = min(ratingscale_max, predict);
        predict = max(ratingscale_min, predict);
    end
end

end

