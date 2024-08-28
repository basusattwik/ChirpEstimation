function updateLivePlots(obj, tind, nind, ax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for pind = 1:obj.numParticles

    objp    = reshape(squeeze(obj.saveObjFunc(pind,1:tind, 1:nind)), tind * nind, 1);
    avgObjp = objp;%smoothdata(objp, 'sgolay', 20);
    set(ax(pind), 'XData', 1:tind*nind, 'YData', avgObjp);
    set(ax(pind), 'DisplayName', ['particle ', num2str(pind)]);
    lgd = legend('show', 'Location', 'best');
    fontsize(lgd, 12, 'points');

end
end