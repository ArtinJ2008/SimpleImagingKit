using Test
using SimpleImagingKit

@testset "windowing" begin
    img = fill(0.0, 4, 4)
    WW, WL = SimpleImagingKit.window_preset(:brain)
    out = SimpleImagingKit.apply_window(img, WW, WL)
    @test size(out) == size(img)
    @test 0.0 <= minimum(out) <= maximum(out) <= 1.0
end
