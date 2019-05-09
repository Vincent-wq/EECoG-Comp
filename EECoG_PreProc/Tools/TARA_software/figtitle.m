function tt = figtitle(str)

% tt = figtitle(str)

ax = axes('position', [0 0 1 1]);
tt = text(0.5, 1.0, str, 'color', 'black', 'units', 'normalized', 'horizontalalignment', 'center', 'verticalalignment', 'top');
set(ax,'visible', 'off')

