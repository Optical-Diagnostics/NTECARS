function spectra_compariso_with_residual(s1::Spectrum, s2::Spectrum; label1="spectrum_1", label2="spectrum_1")
    fig = Figure(size = (600, 500))
    ax = Axis(fig[1,1], yticks = 0:0.1:1)

    lines!(ax, 1e9 .* wavelengths(s1), (intensities(s1, normalization= :maximum)), label = label1, color = :black)
    lines!(ax, 1e9 .* wavelengths(s2), (intensities(s2, normalization= :maximum)), label = label2, color = :red, linestyle = :dash)
    axislegend(ax)
    fig
end