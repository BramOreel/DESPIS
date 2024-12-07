function [bins_a,bins_b] = fixed_transmitter_side_beamformer(h1,h2,N)

%Convert h1 & h2 to the frequency domain
H1 = fft(h1,N/2-1);
H2 = fft(h2,N/2-1);
H12 = sqrt(H1.*conj(H1) + H2.*conj(H2));

bins_a = (conj(H1)./H12);
bins_b = (conj(H2)./H12);

gen_fig = 0;

if(gen_fig)
    figure
    plot(abs(H1(10:end-10)).^2);
    hold on
    plot(abs(H2(10:end-10)).^2);
    plot(abs(H12(10:end-10)).^2);
    legend('H1','H2','H3');
    title('Frequency responses H1, H2, H1+2');
    hold off
end

