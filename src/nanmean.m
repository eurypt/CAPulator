%{
Whe plottin histograms / bar plots, gramm uses nanmean from the statistics
toolbox. Thus, I redefined nanmean here in case the user does not have the
statistics toolbox.
%}
function y = nanmean(x)

y = mean(x(~isnan(x)));