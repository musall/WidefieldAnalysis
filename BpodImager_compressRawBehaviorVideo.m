sPath = 'X:\smusall\BehaviorVideo\Animals\';

animals = dir(sPath);
for iAnimals = 3 : length(animals)
    
    Behavior_compressRawVideo([sPath animals(iAnimals).name filesep], animals(iAnimals).name)
    
end