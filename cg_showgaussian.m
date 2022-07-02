function cg_showgaussian

matname = spm_select(1,'.*_sn.mat','Select parameter file');

sn = load(matname);

Image = spm_read_vols(sn.VF);
[H,X] = hist(Image(:),256);
H = H(2:255);
X = X(2:255);
nH = H/sum(H(:));

nclass = sum(sn.flags.ngaus(1:3));

G = zeros(nclass,length(X));
for i=1:nclass
	G(i,:) = 5*gauss(X,sn.flags.mn(i),sn.flags.vr(i))*sn.flags.mg(i);
end

figure(11)
% Plot the fitted gaussians and boundaries
plot(X,nH,':');
hold on;
for i=1:nclass
  plot(X,G(i,:),'--');
end
plot(X,sum(G),'-');
hold off;


function y=gauss(cr,mn,vr)

y = exp((cr-mn).^2/(-2*vr))/sqrt(2*pi*vr+eps);

return