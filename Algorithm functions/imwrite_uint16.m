function [] = imwrite_uint16(img,fn)
    img = double(img);
    img = img-min(min(img));
    img = img./max(max(img)).*65535;
    img = uint16(img);
    t = Tiff(fn, 'w'); 
    tagstruct.ImageLength = size(img, 1); 
    tagstruct.ImageWidth = size(img, 2); 
    tagstruct.Compression = Tiff.Compression.Deflate; 
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt; 
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack; 
    tagstruct.BitsPerSample = 16; 
    tagstruct.SamplesPerPixel = 1; 
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; 
    t.setTag(tagstruct); 
    t.write(img); 
    t.close();
end