
for npasses = 1:75
    figure(83)
    subplot(1,2,1)
    plot(rad_vs_beam(:,npasses))
    hold on
    plot(profile_b*10e10,'r')
    hold off
    ylim([0 25e10])
    subplot(1,2,2)
    plot(max(rad_vs_und(:,1:npasses)),'b')
    xlim([0 75])
    pause(0.1)
    M(npasses) = getframe(gcf);
end
v = VideoWriter('newfile.avi');
v.FrameRate = 5;
open(v)
writeVideo(v,M)
close(v)