function [Z,flag]=combine_MDG(data_mat,data_dl,data_goes,iok_mat,iok_dl,iok_goes,opt)
% function [Z]=combine_MDG(data_mat,data_dl,data_goes,iok_mat,iok_dl,iok_goes,opt);
% choose option 1 if the datalogger diagnostics are available
% choose option 2 if the are not. Option 2 will use the goes diagnostics for
% the data logger data

if opt==1
    Z = NaN*ones(1,length(data_mat));
    flag = NaN*ones(1,length(data_mat));

    Z (iok_mat) = data_mat(iok_mat);
    flag(iok_mat) = 1;

    Z (~iok_mat & iok_dl) = data_dl(~iok_mat & iok_dl);
    flag(~iok_mat & iok_dl)=2;

    Z (~iok_mat & ~iok_dl & iok_goes) = data_goes(~iok_mat & ~iok_dl & iok_goes);
    flag(~iok_mat & ~iok_dl & iok_goes)=3;
elseif opt==2
    Z = NaN*ones(1,length(data_mat));
    Z (iok_mat) = data_mat(iok_mat);
    Z (~iok_mat & iok_goes) = data_dl(~iok_mat & iok_goes);
    Z (~iok_mat & iok_goes) = data_goes(~iok_mat & iok_goes);
end 
return
    
    


