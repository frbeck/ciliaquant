function [ciliaI] = artifactCleaner(ciliaI0)

grayImage = ciliaI0;
avgBrightness = mean2(grayImage);

mask = grayImage > 3*avgBrightness;
grayImage(~mask) = 0;
ciliaI = grayImage;

end


