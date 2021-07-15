% Function to compute a weighted histogram of 2d data
% data: two dimensional particle data array
% weights: weight of each particle
% nbins: No. of bins in the histogram
% plot: if value is 1 a contour plot of the histogram is displayed
% filter: Binomial filter is applied to the histogram data "filter" times
%
% author: Katharina Kormann, 2019
function [histw, interval1, interval2] = hist3w(data, weights, nbins,plot,filter)

  min1 = min(data(:,1));
  max1 = max(data(:,1));
  delta1 = (max1-min1)/nbins(1);
  interval1 = linspace(min1, max1, nbins(1))-delta1*0.5;
 
  min2 = min(data(:,2));
  max2 = max(data(:,2));
  delta2 = (max2-min2)/nbins(2);
  interval2 = linspace(min2, max2, nbins(2))-delta2*0.5;
  
  histw = zeros(nbins(1),nbins(2));
  
  for i=1:length(data)
      ind(1) = find(interval1 < data(i,1), 1, 'last');
      if ~isempty(ind(1))
          ind(2) = find(interval2 < data(i,2),1, 'last');
          if ~isempty(ind(2))
            histw(ind(1),ind(2)) = histw(ind(1),ind(2)) + weights(i);
          end
      end
  end
  
  histw = histw/(delta1*delta2);
  
  if (filter > 0)
      for j=1:filter
        data2 = histw;
        histw = filter_periodic(data2);
      end
  end
  
  if (plot == 1)
    [xx,vv] = ndgrid(interval1,interval2);
    contour(xx,vv,histw);
    xlabel('x');ylabel('v');
  end
