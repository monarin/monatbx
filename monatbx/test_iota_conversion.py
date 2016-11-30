import prime.iota.iota_conversion as i2p

if (__name__ == "__main__"):
  beam_center = [2052.50, 2078.45]
  beamstop = 0
  square_mode = 'crop'
  img_in = "/reg/d/psdm/xpp/xppg9515/usr/marccd/brunger/397/A5/36-1_0_00161.mccd"
  img_out = "36-1_0_00161.pickle"
  i2p.convert_image(img_in, img_out, square_mode, beamstop, beam_center)
