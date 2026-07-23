using VlasovMaxwellDispersion
using ReTestItems

const NW = parse(Int, get(ENV, "TEST_NWORKERS", "0"))

if NW == 0
    runtests(VlasovMaxwellDispersion; nworkers=0, testitem_timeout=1800)
else
    runtests(ti -> !(:latency in ti.tags), VlasovMaxwellDispersion; nworkers=NW, testitem_timeout=1800)
    runtests(ti -> (:latency in ti.tags), VlasovMaxwellDispersion; nworkers=0, testitem_timeout=1800)
end
