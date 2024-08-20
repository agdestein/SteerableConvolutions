rotmat(angle) = [cos(angle) -sin(angle); sin(angle) cos(angle)]
function flipmat(flip)
    s = flip ? -1 : 1
    [1 0; 0 s]
end
function fliprotmat(flip, angle)
    s = flip ? -1 : 1
    [cos(angle) -s * sin(angle); sin(angle) s * cos(angle)]
end
