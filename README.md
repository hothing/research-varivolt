# VARIX controller

This controller is used for variable voltage transformator TAMINI-TES (I called it "VARIX"). This transformer has big advantage - there is no sliding (mechanical) contact in the secondary windings. 

The main demand is very accurate voltage regulation without windings position encoders 
and with very noised voltage measurement.

After investigation I found that the best controller is state-space controller with Kalman observer (filter).

## Implementations

### Parameters calculation
See a script 'varix-calc-params.sce' is using to get parameters.
The Scilab has been used as engineering enviroment.

Input:
* minimal voltage [V]
* maximal voltage [V]
* minimal position of windings [mm]
* maximal position of windings [mm]
* translation factor [mm/360 deg]
* motor/drive time lag [s]
* nominal motor speed [1/min]

Output:
* dynamic model maitrices : A, B, C
* Kalman observer/filter gain : Lk
* SS-controller feedback gain matrix and setpoint gain: Ks, Nc 

### Siemens TIA Portal

The controller has been implemented for the S7-1500 (TIA-Portal V14).
See file 'VARI_CONTX2.scl'.
