function out = bs(cf,w)
%this function takes a wealth vector and coefficients and returns budget shares

cf_mat = repmat(cf,[1,1,size(w,1)]);
w_mat = repmat([ones(size(w,1),1),w,w.^2]',[1,1,28]);
w_mat = permute(w_mat,[3 1 2]);
log_list = squeeze(log(sum(cf_mat.*w_mat,2)));
out_part = bsxfun(@rdivide,log_list',sum(log_list',2)+1);

out = [1-sum(out_part,2),out_part];

end