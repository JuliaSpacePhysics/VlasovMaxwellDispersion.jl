module VlasovMaxwellDispersionUnitfulExt

import VlasovMaxwellDispersion
import VlasovMaxwellDispersion.PlasmaBase as PB

include(joinpath(pkgdir(VlasovMaxwellDispersion), "lib", "PlasmaBase", "src", "unitful.jl"))

end
