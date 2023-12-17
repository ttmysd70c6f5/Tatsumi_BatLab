function flight_segmentation(sessionpath)
load(sessionpath)

%=== FIGURE: Velocity and flight segmentation
figure(); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
tiledlayout(n_tags,1,'TileSpacing','tight');
for i = 1:n_tags
    ax(i) = nexttile;   %ax(i) = subplot(n_tags,1,i);
    area(t,bflying(:,i)*5,'FaceAlpha',0.3,'LineStyle','none');  hold on;
    area(t,wBeats(:,i)*-1,'FaceAlpha',0.3,'LineStyle','none');  refline(0,-0.3);
    plot(t,v_abs(:,i),'.','Color', bat_clr(i,:));     plot(t,r(:,1,i),'k--');  ylabel('Velocity (m/s)');     hold off;
    legend('Fly','Wing-B','Vel','x(m)');
    title([num2str(f_num(i)) ' flights']);
end
linkaxes(ax,'x');   xlabel('Time(s)');

end