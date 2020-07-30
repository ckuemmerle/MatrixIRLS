function predict = get_prediction_RecSys_IRLS(M,rowind_test, ...
    colind_test,centscale,clip_ratingscale,training_flag)

m=length(rowind_test);
data=struct;
data.rowind = rowind_test;
data.colind = colind_test;
data.d1     = size(M.U,1);
data.d2     = size(M.V,1);

predict = partXY((M.U*M.Gam1+M.Gam3)',M.V',rowind_test,colind_test,m)';
predict = predict+partXY(M.U.',M.Gam2,rowind_test,colind_test,m)';
if training_flag
    predict = predict+M.res_range.';
end

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

