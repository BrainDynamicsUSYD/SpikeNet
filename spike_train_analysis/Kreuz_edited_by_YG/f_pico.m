function [ave,plot_x_values,plot_y_values]=f_pico(isi,pico,dt,tmin)

cum_isi=cumsum([0 isi]);
ave=sum(pico.*repmat(isi,size(pico,1),1),2)/cum_isi(end);
plot_x_values=tmin+unique([cum_isi(1:end) cum_isi(2:end-1)-dt/100]);
plot_y_values=reshape([pico; pico],1,length(pico)*2);

