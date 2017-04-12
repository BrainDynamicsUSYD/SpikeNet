function R = get_dist_CC(R, n_p)
% distance-dependent spike count correlation coefficient
pair = R.Analysis.CC_pop_pairs{1};
cc = R.Analysis.CC_pop{1};

if nargin < 2
    n_p = 20;
end

fw = sqrt(R.N(1));
hw = (fw-1)/2;
if mod(hw,1) ~= 0
    error('Not a supported grid.')
end
Lat = lattice_nD(2,hw);
Lat = Lat';



p_A = Lat(:,pair(1,:));
p_B = Lat(:,pair(2,:));
dX = min(abs(p_A(1,:)-p_B(1,:)), fw - abs(p_A(1,:)-p_B(1,:)));
dY = min(abs(p_A(2,:)-p_B(2,:)), fw - abs(p_A(2,:)-p_B(2,:)));
d = sqrt( dX.^2 + dY.^2 );

d_bin = linspace(0, max(d), n_p);
[meanData, bin_mid] = histcn(d(:), d_bin,  'AccumData', cc(:), 'Fun', @nanmean);
stdData = histcn(d(:), d_bin,  'AccumData', cc(:), 'Fun', @nanstd);


R.Analysis.dist_cc.d_bin = d_bin;
R.Analysis.dist_cc.mu = meanData;
R.Analysis.dist_cc.sig = stdData;

end
