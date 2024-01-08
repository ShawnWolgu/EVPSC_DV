# Relative equations
## Definitions
The polycrystal thermal-elasticity strain

$$ \boldsymbol E_t = \boldsymbol\alpha \delta T$$ 

The temperature diffence with respect the reference temperature

$$ \delta T(t) = T(t) - T_{ref}(t) $$

The reference temperature is the temperature with free of thermal stress.

## In thermal-elasticity self-consistency
$$ \boldsymbol B^e=\left(\boldsymbol M^e + \widetilde{\boldsymbol M}^e\right)^{-1}:\left(\overline{\boldsymbol M}^e+\widetilde{\boldsymbol M}^e\right) $$

$$ \overline{\boldsymbol M}^e=\langle\boldsymbol M^e:\boldsymbol B^e\rangle\:\langle\boldsymbol B^e\rangle^{-1} $$

$$ \boldsymbol B^{e+}=\left(\overline{\boldsymbol M}^e+\widetilde{\boldsymbol M}^e\right):\left(\boldsymbol M^e + \widetilde{\boldsymbol M}^e\right)^{-1} $$

$$ \overline{\boldsymbol \alpha}=\langle\boldsymbol B^{e+}\rangle^{-1}:\langle\boldsymbol B^{e+} :\boldsymbol \alpha\rangle\ $$

## In single crystal

$$ \widetilde{\boldsymbol d} = -\widetilde{\boldsymbol M}^e:(\dot{\boldsymbol\sigma^\nabla}-\dot{\boldsymbol\Sigma^\nabla})-\widetilde{\boldsymbol M}^{vp}:(\boldsymbol\sigma-\boldsymbol\Sigma) $$

$$ \widetilde{\boldsymbol d}=\boldsymbol d-\overline{\boldsymbol d}=\boldsymbol d_e + \boldsymbol d_p+\boldsymbol d_t -\overline{\boldsymbol d}$$

$$ \boldsymbol d_t=\boldsymbol\alpha\frac{d \delta T}{dt} = \boldsymbol\alpha\frac{dT}{dt} $$

## In polycrystal

$$ \boldsymbol D =\overline{\boldsymbol M}^e:\dot{\boldsymbol\Sigma^\nabla}+\boldsymbol D^t+\overline{\boldsymbol M}^{vp}:\boldsymbol\Sigma+\boldsymbol D^0=\langle\boldsymbol M^e:\boldsymbol B^e\rangle:\boldsymbol\Sigma^\nabla+\dot{\boldsymbol E_t}+\langle\boldsymbol M^{vp}:\boldsymbol B^{vp}\rangle:\boldsymbol\Sigma+\langle\boldsymbol M^{vp}:\boldsymbol b^{vp}+\boldsymbol d^0\rangle $$
