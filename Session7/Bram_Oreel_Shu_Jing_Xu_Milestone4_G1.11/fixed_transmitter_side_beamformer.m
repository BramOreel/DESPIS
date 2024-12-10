function [bins_a,bins_b] = fixed_transmitter_side_beamformer(H)

%H has to be size N/2-1
H1 = H(:,1);
H2 = H(:,2);
H1_full = [0;H1 ;0; flipud(conj(H1))];
H2_full = [0;H2 ;0; flipud(conj(H2))];
H12 = sqrt(H1.*conj(H1) + H2.*conj(H2));
H12_full = [0;H12 ;0; flipud(conj(H12))];
bins_a = (conj(H1)./H12);
bins_b = (conj(H2)./H12);

gen_fig = 1;

if(gen_fig)
    figure
    plot(abs(H1_full(10:end-10)).^2);
    hold on
    plot(abs(H2_full(10:end-10)).^2);
    plot(abs(H12_full(10:end-10)).^2);
    legend('H1','H2','H3');
    title('Frequency responses H1, H2, H1+2');
    hold off
end

