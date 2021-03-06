function [] = imwrite_float(img,fn)
    t = Tiff(fn, 'w'); 
    tagstruct.ImageLength = size(img, 1); 
    tagstruct.ImageWidth = size(img, 2); 
    tagstruct.Compression = Tiff.Compression.Deflate; 
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP; 
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
    tagstruct.BitsPerSample = 32; 
    tagstruct.SamplesPerPixel = 1; 
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
    t.setTag(tagstruct); 
    t.write(img); 
    t.close();
end