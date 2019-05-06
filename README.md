# pcamod
Data-driven process modeling for monitoring purposes using PCA

Ejemplo:

```[MATLAB]
% Entrenamiento fuera de línea
X = [...];                        % Datos de entrenamiento
model = pcamod(X);                % Modelo del comportamiento normal 
[uT2,uSPE] = model.ucl(3,0.95);   % Hallar umbrales al 95% con 3 PCs
% Operación en línea
xnew = [...];                     % Dato de prueba
[~,~,T2,SPE] = model.map(xnew,3); % Calcular los índices T^2 y SPE
if T2 > uT2 | SPE > uSPE          % ¿Algún índice excede el umbral?
   disp('¡Falla!')                % ...se infiere que existe falla.
end
```
