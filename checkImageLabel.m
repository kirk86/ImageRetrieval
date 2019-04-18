function label = checkImageLabel(imageName)
% return the equivalent label of each image

if (imageName >= 0 && imageName <= 99)
    label = 1;
elseif (imageName > 99 && imageName <= 199)
    label = 2;
elseif (imageName > 199 && imageName <= 299)
    label = 3;
elseif (imageName > 299 && imageName <= 399)
    label = 4;
elseif (imageName > 399 && imageName <= 499)
    label = 5;
elseif (imageName > 499 && imageName <= 599)
    label = 6;
elseif (imageName > 599 && imageName <= 699)
    label = 7;
elseif (imageName > 699 && imageName <= 799)
    label = 8;
elseif (imageName > 799 && imageName <= 899)
    label = 9;
elseif (imageName > 899 && imageName <= 999)
    label = 10;
end

end