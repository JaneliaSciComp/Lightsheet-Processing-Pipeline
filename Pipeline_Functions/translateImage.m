function im_out = translateImage( im_in, trans_x, trans_y, trans_z )
im_out = imtranslate(im_in, [trans_x, trans_y, trans_z]);
end

