clear all; close all;

addpath('events');

typeA = eventType.event

typeB = eventType.collision

if typeA ~= typeB
    display('NOPE');
end
