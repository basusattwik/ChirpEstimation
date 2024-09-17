function updateLivePlots_SL(obj, tind, ax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for pind = 1:obj.numParticles

    idx  = 1:tind;
    objp = obj.saveObjFunc(pind, idx);
    set(ax(pind), 'XData', idx, 'YData', objp);
    set(ax(pind), 'DisplayName', ['particle ', num2str(pind)]);
    lgd = legend('show', 'Location', 'best');
    fontsize(lgd, 12, 'points');

end

end