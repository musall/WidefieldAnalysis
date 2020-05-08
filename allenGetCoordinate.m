function allenGetCoordinate(dorsalMaps)
% code to find anatomical coordinates in the allen framework.
% reports coordinates in mm, relative to bregma location.

if ~exist('dorsalMaps','var') || isempty(dorsalMaps)
    load allenDorsalMapSM.mat
end

figure('Renderer','painters');
cImg = imagesc(dorsalMaps.dorsalMapScaled,'ButtonDownFcn',{@showCoords,dorsalMaps}); axis image; colormap colorcube
hold on; 
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{x}(:,2), dorsalMaps.edgeOutlineSplit{x}(:,1), 'w', 'LineWidth', 0.2,'ButtonDownFcn',{@showCoords,dorsalMaps});axis image
end
plot(dorsalMaps.bregmaScaled(1),dorsalMaps.bregmaScaled(2),'k.','MarkerSize',20,'linewidth',2,'ButtonDownFcn',{@showCoords,dorsalMaps}); %bregma location
ylim([1 550])

pxPerMM = dorsalMaps.pxPerMM;
ax = cImg.Parent;

ax.XTick = [fliplr(dorsalMaps.bregmaScaled(1) : -dorsalMaps.pxPerMM : 0) dorsalMaps.bregmaScaled(1) + pxPerMM : pxPerMM : size(dorsalMaps.dorsalMapScaled,2)];
ax.XTickLabels = -floor(dorsalMaps.bregmaScaled(1) / pxPerMM) :  length(ax.XTick) - floor(dorsalMaps.bregmaScaled(1) / pxPerMM);

ax.YTick = [fliplr(dorsalMaps.bregmaScaled(2) : -dorsalMaps.pxPerMM : 0) dorsalMaps.bregmaScaled(2) + pxPerMM : pxPerMM : size(dorsalMaps.dorsalMapScaled,1)];
ax.YTickLabels = -floor(dorsalMaps.bregmaScaled(2) / pxPerMM) :  length(ax.XTick) - floor(dorsalMaps.bregmaScaled(2) / pxPerMM);
grid on


end

function showCoords(hObject,~,dorsalMaps)
   ax = hObject.Parent;
   pointerSize = 5;
   cCoords = get(gca,'CurrentPoint');
   cVal = dorsalMaps.dorsalMapScaled(round(cCoords(1,2)), round(cCoords(1,1)));
   if cVal > 0
       cArea = dorsalMaps.labelTable.name{find(dorsalMaps.labelTable.id == cVal, 1)};
       cCode = dorsalMaps.labelTable.abbreviation{find(dorsalMaps.labelTable.id == cVal, 1)};
       delete(findall(ax.Children,'Type','rectangle'));
       Frame = [cCoords(1,1)-pointerSize cCoords(1,2)-pointerSize pointerSize*2 pointerSize*2];
       rectangle('Position',Frame,'linewidth',2,'edgecolor','k','Curvature',1,'parent',ax,'ButtonDownFcn',{@showCoords,dorsalMaps});
       cCoords = round((cCoords(1,1:2) - dorsalMaps.bregmaScaled) ./ dorsalMaps.pxPerMM, 2);
       fprintf('Current area: %s (%s)\n%g mm ML; %g mm RC\n', cArea, cCode, abs(cCoords(1)), -cCoords(2));
   end
end