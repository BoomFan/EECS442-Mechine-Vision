function [points1 points2] = readPoints(file)
    points = textread(file);
    points1 = points(:, 1:2);
    points2 = points(:, 3:4);
end