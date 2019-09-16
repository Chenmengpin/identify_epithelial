temp_png_function <-
  function(x) {
    png(
      file = (x), 
      width = 14, 
      height = 8, 
      res = 300, 
      units = 'in'
    )
  }
