using Test
import SimpleImagingKit.Series as SIK

@testset "windowing" begin
    img = fill(0.0, 4, 4)
    WW, WL = SIK.window_preset(:brain)
    out = SIK.apply_window(img, WW, WL)
    @test size(out) == size(img)
    @test 0.0 <= minimum(out) <= maximum(out) <= 1.0
end
