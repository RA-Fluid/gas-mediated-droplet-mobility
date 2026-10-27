function blendedRGB = RA_getAlphaColor(color, alpha, bgColor)
    % color: 'hex' string or [R G B] vector (0-1)
    % alpha: transparency value from 0 to 1
    % bgColor: [R G B] of your plot background (default is white [1 1 1])

    if nargin < 3, bgColor = [1 1 1]; end
    
    % Convert Hex to RGB if necessary
    if ischar(color) || isstring(color)
        if startsWith(color, '#'), color = erase(color, '#'); end
        fgRGB = reshape(sscanf(color, '%2x'), 1, 3) / 255;
    else
        fgRGB = color;
    end

    % Ensure alpha is within [0, 1]
    alpha = max(0, min(1, alpha));

    % Apply Alpha Blending Formula
    blendedRGB = alpha * fgRGB + (1 - alpha) * bgColor;
end