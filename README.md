# VARIX controller

This controller is used for variable voltage transformator TAMINI-TES (I called it "VARIX"). This transformer has big advantage - there is no sliding (mechanical) contact in the secondary windings. 

The main demand is very accurate voltage regulation without windings position encoders 
and with very noised voltage measurement.

After investigation I found that the best controller is state-space controller with Kalman observer (filter).
