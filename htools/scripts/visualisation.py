import yt

def vol_render(de, fname, res=1024, bmin = 1e-6,bmax = 1e3):
    data = dict(density = (de, "g/cm**3"))
    ds = yt.load_uniform_grid(data, de.shape, length_unit="pc")
    sc = yt.create_scene(ds, field=("density"))
    sc.camera.resolution = (res, res)
    sc.camera.focus = ds.arr([0.3, 0.3, 0.3], "unitary")
    source = sc[0]
    source.tfh.set_bounds((bmin, bmax))
    sc.camera.position = ds.arr([0, 0, 0], "unitary")
    sc.render()
    sc.save(sim_dir+'/'+fname+'.png', sigma_clip=4)
