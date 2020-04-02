# 物理常数

This chapter describes macros for the values of physical constants, such as the speed of light, ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png), and gravitational constant, ![G](https://www.gnu.org/software/gsl/doc/html/_images/math/a205209d309440659169ad888f53d45f717d8d0a.png). The values are available in different unit systems, including the standard MKSA system (meters, kilograms, seconds, amperes) and the CGSM system (centimeters, grams, seconds, gauss), which is commonly used in Astronomy.

本章描述物理常数的值的宏，比如光速$$c$$，和引力常数，$$G$$。可以获得不同单位体系下的这些值，包括标准MKSA体系(米，千克，秒，安倍)和CGSM体系(厘米，克，秒，高斯)，后者通常在天文学中用到。

The definitions of constants in the MKSA system are available in the file `gsl_const_mksa.h`. The constants in the CGSM system are defined in `gsl_const_cgsm.h`. Dimensionless constants, such as the fine structure constant, which are pure numbers are defined in `gsl_const_num.h`.

在MKSA体系下的常数的定义在文件`gsl_const_mksa.h`中可以找到。在CGSM体系下的常数定义在文件`gsl_const_cgsm.h`中。无维的常数，比如精细结构常数，通常为纯数字被定义在`gsl_const_num.h`中。

The full list of constants is described briefly below. Consult the header files themselves for the values of the constants used in the library.

常数的完整列表在下面得到了简短描述。查询头文件本身以知晓库中使用的常数的值。



## 基本常数

- `GSL_CONST_MKSA_SPEED_OF_LIGHT`

  The speed of light in vacuum, ![c](https://www.gnu.org/software/gsl/doc/html/_images/math/a2ff0b909e008fecae54f397fc3bd6ef4ae96b5c.png).

  真空中的光速，$$c$$。

- `GSL_CONST_MKSA_VACUUM_PERMEABILITY`

  The permeability of free space, ![\mu_0](https://www.gnu.org/software/gsl/doc/html/_images/math/5aa6e1c125e18e7278a43635e3f703e590567f7a.png). This constant is defined in the MKSA system only.

  自由空间的磁导率，$$\mu_0$$。这个常数仅以MKSA体系定义。

- `GSL_CONST_MKSA_VACUUM_PERMITTIVITY`

  The permittivity of free space, ![\epsilon_0](https://www.gnu.org/software/gsl/doc/html/_images/math/ff8259230190e6a23c5c3db6c2874cf092dfa597.png). This constant is defined in the MKSA system only.

  自由空间的介电常数，$$\epsilon_0$$。这个常数仅以MKSA体系定义。

- `GSL_CONST_MKSA_PLANCKS_CONSTANT_H`

  Planck’s constant, ![h](https://www.gnu.org/software/gsl/doc/html/_images/math/40537c284ff40128b45f39365983e63b41d9f731.png).

  普朗克常数，$$h$$。

- `GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR`

  Planck’s constant divided by ![2\pi](https://www.gnu.org/software/gsl/doc/html/_images/math/059fbd5d7d78643be5cf4a37b45b0f8140e2d605.png), ![\hbar](https://www.gnu.org/software/gsl/doc/html/_images/math/744f2e3772d21829d5709397eb94e48454abe400.png).

  普朗克常数除以$$2\pi$$，$$\hbar$$。

- `GSL_CONST_NUM_AVOGADRO`

  Avogadro’s number, ![N_a](https://www.gnu.org/software/gsl/doc/html/_images/math/43e1b83f4660f134450bf96ea213546bd11710fc.png).

  阿伏伽德罗常数，$$N_a$$。

- `GSL_CONST_MKSA_FARADAY`

  The molar charge of 1 Faraday.

  1法拉第的摩尔电荷。

- `GSL_CONST_MKSA_BOLTZMANN`

  The Boltzmann constant, ![k](https://www.gnu.org/software/gsl/doc/html/_images/math/90645b5c2abc13d7de104fd87abeb2f19406e5f8.png).

  玻尔兹曼常数，$$k$$。

- `GSL_CONST_MKSA_MOLAR_GAS`

  The molar gas constant, ![R_0](https://www.gnu.org/software/gsl/doc/html/_images/math/b2ca85120822033eb09f9a1b72245dd537afdbaa.png).

  摩尔气体常数，$$R_0$$。

- `GSL_CONST_MKSA_STANDARD_GAS_VOLUME`

  The standard gas volume, ![V_0](https://www.gnu.org/software/gsl/doc/html/_images/math/00136b16b615351d0f38ce9132fae0ed3b00fa61.png).

  标准气体体积，$$V_0$$。

- `GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT`

  The Stefan-Boltzmann radiation constant, ![\sigma](https://www.gnu.org/software/gsl/doc/html/_images/math/d4e96af4d650f28ddd3f3959a0c1054d4098ed8f.png).

  斯蒂芬-玻尔兹曼辐射常数，$$\sigma$$。

- `GSL_CONST_MKSA_GAUSS`

  The magnetic field of 1 Gauss.

  1高斯的磁场。



## Astronomy and Astrophysics

- `GSL_CONST_MKSA_ASTRONOMICAL_UNIT`

  The length of 1 astronomical unit (mean earth-sun distance), ![au](https://www.gnu.org/software/gsl/doc/html/_images/math/94238ed6e0d720f92f2ce5270f546ae82a666e8b.png).

- `GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT`

  The gravitational constant, ![G](https://www.gnu.org/software/gsl/doc/html/_images/math/a205209d309440659169ad888f53d45f717d8d0a.png).

- `GSL_CONST_MKSA_LIGHT_YEAR`

  The distance of 1 light-year, ![ly](https://www.gnu.org/software/gsl/doc/html/_images/math/f033a0f20580c4a9d60fb1709c007fa43c27ffe0.png).

- `GSL_CONST_MKSA_PARSEC`

  The distance of 1 parsec, ![pc](https://www.gnu.org/software/gsl/doc/html/_images/math/5e2678cd6801bda1fe5e087d816af4a04fffad44.png).

- `GSL_CONST_MKSA_GRAV_ACCEL`

  The standard gravitational acceleration on Earth, ![g](https://www.gnu.org/software/gsl/doc/html/_images/math/1954347bc3687d8e02da2d3c02f34acab120853f.png).

- `GSL_CONST_MKSA_SOLAR_MASS`

  The mass of the Sun.



## Atomic and Nuclear Physics

- `GSL_CONST_MKSA_ELECTRON_CHARGE`

  The charge of the electron, ![e](https://www.gnu.org/software/gsl/doc/html/_images/math/9a7ff5d0469d9a20a6a7da1b2c1ab9b1014d53c9.png).

- `GSL_CONST_MKSA_ELECTRON_VOLT`

  The energy of 1 electron volt, ![eV](https://www.gnu.org/software/gsl/doc/html/_images/math/0a0168b8d7c56c4513b58ef24b684486985f6aa7.png).

- `GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS`

  The unified atomic mass, ![amu](https://www.gnu.org/software/gsl/doc/html/_images/math/b551520a748ea76ec63b676668af00a7a3c19b28.png).

- `GSL_CONST_MKSA_MASS_ELECTRON`

  The mass of the electron, ![m_e](https://www.gnu.org/software/gsl/doc/html/_images/math/bb52c6750d0e25caa396c9e745f105e602f7dde4.png).

- `GSL_CONST_MKSA_MASS_MUON`

  The mass of the muon, ![m_\mu](https://www.gnu.org/software/gsl/doc/html/_images/math/e53987ae28004a8313c43263a109ed3f636396a7.png).

- `GSL_CONST_MKSA_MASS_PROTON`

  The mass of the proton, ![m_p](https://www.gnu.org/software/gsl/doc/html/_images/math/788ed659b501fa07f8a956a52671db9974ebb387.png).

- `GSL_CONST_MKSA_MASS_NEUTRON`

  The mass of the neutron, ![m_n](https://www.gnu.org/software/gsl/doc/html/_images/math/ec50182f12d5f365ca0b837fa2691aebf75cc0d6.png).

- `GSL_CONST_NUM_FINE_STRUCTURE`

  The electromagnetic fine structure constant ![\alpha](https://www.gnu.org/software/gsl/doc/html/_images/math/914a439c7f061a4faeb93165a8d4bf26c6688e1a.png).

- `GSL_CONST_MKSA_RYDBERG`

  The Rydberg constant, ![Ry](https://www.gnu.org/software/gsl/doc/html/_images/math/28386cdc38aeb076fca2a398daafeac6fda19056.png), in units of energy. This is related to the Rydberg inverse wavelength ![R_\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/8750ab7132c8b0fd6c6f7398b67fa998fa2d1cb6.png) by ![Ry = h c R_\infty](https://www.gnu.org/software/gsl/doc/html/_images/math/180309f33e20d3f1dfc40b631a3fb93a8b5e202d.png).

- `GSL_CONST_MKSA_BOHR_RADIUS`

  The Bohr radius, ![a_0](https://www.gnu.org/software/gsl/doc/html/_images/math/346baeeeeff5efc6717502bc129c937c17041a0c.png).

- `GSL_CONST_MKSA_ANGSTROM`

  The length of 1 angstrom.

- `GSL_CONST_MKSA_BARN`

  The area of 1 barn.

- `GSL_CONST_MKSA_BOHR_MAGNETON`

  The Bohr Magneton, ![\mu_B](https://www.gnu.org/software/gsl/doc/html/_images/math/e48b87598e9d25e8dedc8c69a5abde620e3af30d.png).

- `GSL_CONST_MKSA_NUCLEAR_MAGNETON`

  The Nuclear Magneton, ![\mu_N](https://www.gnu.org/software/gsl/doc/html/_images/math/3d6a1774bfdf354871ceb78b7eb6e8dc5a00bf37.png).

- `GSL_CONST_MKSA_ELECTRON_MAGNETIC_MOMENT`

  The absolute value of the magnetic moment of the electron, ![\mu_e](https://www.gnu.org/software/gsl/doc/html/_images/math/c0b70d281450571c683f1a2d4f13aa867e8b0066.png). The physical magnetic moment of the electron is negative.

- `GSL_CONST_MKSA_PROTON_MAGNETIC_MOMENT`

  The magnetic moment of the proton, ![\mu_p](https://www.gnu.org/software/gsl/doc/html/_images/math/14a12bec7e880b31c25aebd0816e1c1da67f1a22.png).

- `GSL_CONST_MKSA_THOMSON_CROSS_SECTION`

  The Thomson cross section, ![\sigma_T](https://www.gnu.org/software/gsl/doc/html/_images/math/d881e0d0bc2f1cdac3ab95244b282dd4bd3c813d.png).

- `GSL_CONST_MKSA_DEBYE`

  The electric dipole moment of 1 Debye, ![D](https://www.gnu.org/software/gsl/doc/html/_images/math/0638e3c183dd8cfd1790c25e762f1a8c53ee5711.png).



## Measurement of Time

- `GSL_CONST_MKSA_MINUTE`

  The number of seconds in 1 minute.

- `GSL_CONST_MKSA_HOUR`

  The number of seconds in 1 hour.

- `GSL_CONST_MKSA_DAY`

  The number of seconds in 1 day.

- `GSL_CONST_MKSA_WEEK`

  The number of seconds in 1 week.



## Imperial Units

- `GSL_CONST_MKSA_INCH`

  The length of 1 inch.

- `GSL_CONST_MKSA_FOOT`

  The length of 1 foot.

- `GSL_CONST_MKSA_YARD`

  The length of 1 yard.

- `GSL_CONST_MKSA_MILE`

  The length of 1 mile.

- `GSL_CONST_MKSA_MIL`

  The length of 1 mil (1/1000th of an inch).



## Speed and Nautical Units

- `GSL_CONST_MKSA_KILOMETERS_PER_HOUR`

  The speed of 1 kilometer per hour.

- `GSL_CONST_MKSA_MILES_PER_HOUR`

  The speed of 1 mile per hour.

- `GSL_CONST_MKSA_NAUTICAL_MILE`

  The length of 1 nautical mile.

- `GSL_CONST_MKSA_FATHOM`

  The length of 1 fathom.

- `GSL_CONST_MKSA_KNOT`

  The speed of 1 knot.



## Printers Units

- `GSL_CONST_MKSA_POINT`

  The length of 1 printer’s point (1/72 inch).

- `GSL_CONST_MKSA_TEXPOINT`

  The length of 1 TeX point (1/72.27 inch).



## Volume, Area and Length

- `GSL_CONST_MKSA_MICRON`

  The length of 1 micron.

- `GSL_CONST_MKSA_HECTARE`

  The area of 1 hectare.

- `GSL_CONST_MKSA_ACRE`

  The area of 1 acre.

- `GSL_CONST_MKSA_LITER`

  The volume of 1 liter.

- `GSL_CONST_MKSA_US_GALLON`

  The volume of 1 US gallon.

- `GSL_CONST_MKSA_CANADIAN_GALLON`

  The volume of 1 Canadian gallon.

- `GSL_CONST_MKSA_UK_GALLON`

  The volume of 1 UK gallon.

- `GSL_CONST_MKSA_QUART`

  The volume of 1 quart.

- `GSL_CONST_MKSA_PINT`

  The volume of 1 pint.



## Mass and Weight

- `GSL_CONST_MKSA_POUND_MASS`

  The mass of 1 pound.

- `GSL_CONST_MKSA_OUNCE_MASS`

  The mass of 1 ounce.

- `GSL_CONST_MKSA_TON`

  The mass of 1 ton.

- `GSL_CONST_MKSA_METRIC_TON`

  The mass of 1 metric ton (1000 kg).

- `GSL_CONST_MKSA_UK_TON`

  The mass of 1 UK ton.

- `GSL_CONST_MKSA_TROY_OUNCE`

  The mass of 1 troy ounce.

- `GSL_CONST_MKSA_CARAT`

  The mass of 1 carat.

- `GSL_CONST_MKSA_GRAM_FORCE`

  The force of 1 gram weight.

- `GSL_CONST_MKSA_POUND_FORCE`

  The force of 1 pound weight.

- `GSL_CONST_MKSA_KILOPOUND_FORCE`

  The force of 1 kilopound weight.

- `GSL_CONST_MKSA_POUNDAL`

  The force of 1 poundal.



## Thermal Energy and Power

- `GSL_CONST_MKSA_CALORIE`

  The energy of 1 calorie.

- `GSL_CONST_MKSA_BTU`

  The energy of 1 British Thermal Unit, ![btu](https://www.gnu.org/software/gsl/doc/html/_images/math/02c57b776d130b5e69ec04c3a67e1776b76459ae.png).

- `GSL_CONST_MKSA_THERM`

  The energy of 1 Therm.

- `GSL_CONST_MKSA_HORSEPOWER`

  The power of 1 horsepower.



## Pressure

- `GSL_CONST_MKSA_BAR`

  The pressure of 1 bar.

- `GSL_CONST_MKSA_STD_ATMOSPHERE`

  The pressure of 1 standard atmosphere.

- `GSL_CONST_MKSA_TORR`

  The pressure of 1 torr.

- `GSL_CONST_MKSA_METER_OF_MERCURY`

  The pressure of 1 meter of mercury.

- `GSL_CONST_MKSA_INCH_OF_MERCURY`

  The pressure of 1 inch of mercury.

- `GSL_CONST_MKSA_INCH_OF_WATER`

  The pressure of 1 inch of water.

- `GSL_CONST_MKSA_PSI`

  The pressure of 1 pound per square inch.



## Viscosity

- `GSL_CONST_MKSA_POISE`

  The dynamic viscosity of 1 poise.

- `GSL_CONST_MKSA_STOKES`

  The kinematic viscosity of 1 stokes.



## Light and Illumination

- `GSL_CONST_MKSA_STILB`

  The luminance of 1 stilb.

- `GSL_CONST_MKSA_LUMEN`

  The luminous flux of 1 lumen.

- `GSL_CONST_MKSA_LUX`

  The illuminance of 1 lux.

- `GSL_CONST_MKSA_PHOT`

  The illuminance of 1 phot.

- `GSL_CONST_MKSA_FOOTCANDLE`

  The illuminance of 1 footcandle.

- `GSL_CONST_MKSA_LAMBERT`

  The luminance of 1 lambert.

- `GSL_CONST_MKSA_FOOTLAMBERT`

  The luminance of 1 footlambert.



## Radioactivity

- `GSL_CONST_MKSA_CURIE`

  The activity of 1 curie.

- `GSL_CONST_MKSA_ROENTGEN`

  The exposure of 1 roentgen.

- `GSL_CONST_MKSA_RAD`

  The absorbed dose of 1 rad.



## Force and Energy

- `GSL_CONST_MKSA_NEWTON`

  The SI unit of force, 1 Newton.

- `GSL_CONST_MKSA_DYNE`

  The force of 1 Dyne = ![10^{-5}](https://www.gnu.org/software/gsl/doc/html/_images/math/2ff1650dfefd5c64f9ff5dbd5add4288ff5ec9de.png) Newton.

- `GSL_CONST_MKSA_JOULE`

  The SI unit of energy, 1 Joule.

- `GSL_CONST_MKSA_ERG`

  The energy 1 erg = ![10^{-7}](https://www.gnu.org/software/gsl/doc/html/_images/math/6c6f9131cb760e7552278fb3ce3d0adf0c44d6a4.png) Joule.



## Prefixes

These constants are dimensionless scaling factors.

- `GSL_CONST_NUM_YOTTA`

  ![10^{24}](https://www.gnu.org/software/gsl/doc/html/_images/math/a4438ebf5b77282eab5d2146d68cc757ccfc57e3.png)

- `GSL_CONST_NUM_ZETTA`

  ![10^{21}](https://www.gnu.org/software/gsl/doc/html/_images/math/e00af836c770f1e1fe1de6a7bdc4fc937ede23c2.png)

- `GSL_CONST_NUM_EXA`

  ![10^{18}](https://www.gnu.org/software/gsl/doc/html/_images/math/c54f2c1853d3abe1c28b22cfa113e5728b2d4cde.png)

- `GSL_CONST_NUM_PETA`

  ![10^{15}](https://www.gnu.org/software/gsl/doc/html/_images/math/ad5854921adf6c6ecb6deaea0a593f47a5ca3e3e.png)

- `GSL_CONST_NUM_TERA`

  ![10^{12}](https://www.gnu.org/software/gsl/doc/html/_images/math/e7da6e22cefa892ee0ae5e3d96c6f7baafce12d0.png)

- `GSL_CONST_NUM_GIGA`

  ![10^9](https://www.gnu.org/software/gsl/doc/html/_images/math/c6fd0ac748ce229a171997f88cf5b60cc1d880a0.png)

- `GSL_CONST_NUM_MEGA`

  ![10^6](https://www.gnu.org/software/gsl/doc/html/_images/math/0943135c3f942f9f495c0ebe693a8d02be7f81be.png)

- `GSL_CONST_NUM_KILO`

  ![10^3](https://www.gnu.org/software/gsl/doc/html/_images/math/fdf9cbd4b90ef6b9eca1346cf9c59efe95c53154.png)

- `GSL_CONST_NUM_MILLI`

  ![10^{-3}](https://www.gnu.org/software/gsl/doc/html/_images/math/d10b2ef6892a496823d70981401ea6f6da4bc4a1.png)

- `GSL_CONST_NUM_MICRO`

  ![10^{-6}](https://www.gnu.org/software/gsl/doc/html/_images/math/32189bf572cc2c9d57e219a85372cabe176c2ea4.png)

- `GSL_CONST_NUM_NANO`

  ![10^{-9}](https://www.gnu.org/software/gsl/doc/html/_images/math/b3efae70fa45df5638c5fdeed451b0d8b7393afb.png)

- `GSL_CONST_NUM_PICO`

  ![10^{-12}](https://www.gnu.org/software/gsl/doc/html/_images/math/caa042c6bf911059c20c782df7332a7fb65ea8d6.png)

- `GSL_CONST_NUM_FEMTO`

  ![10^{-15}](https://www.gnu.org/software/gsl/doc/html/_images/math/6c2a0f4b61078df03693839279e17994d3afda34.png)

- `GSL_CONST_NUM_ATTO`

  ![10^{-18}](https://www.gnu.org/software/gsl/doc/html/_images/math/31461163de95c57511f907ef670b8e9b753582f0.png)

- `GSL_CONST_NUM_ZEPTO`

  ![10^{-21}](https://www.gnu.org/software/gsl/doc/html/_images/math/862f268a0d771e872677c128884b3ea39a65bf03.png)

- `GSL_CONST_NUM_YOCTO`

  ![10^{-24}](https://www.gnu.org/software/gsl/doc/html/_images/math/3847094bc15c7f5734bf754148d9593bc13b6428.png)

## Examples

The following program demonstrates the use of the physical constants in a calculation. In this case, the goal is to calculate the range of light-travel times from Earth to Mars.

The required data is the average distance of each planet from the Sun in astronomical units (the eccentricities and inclinations of the orbits will be neglected for the purposes of this calculation). The average radius of the orbit of Mars is 1.52 astronomical units, and for the orbit of Earth it is 1 astronomical unit (by definition). These values are combined with the MKSA values of the constants for the speed of light and the length of an astronomical unit to produce a result for the shortest and longest light-travel times in seconds. The figures are converted into minutes before being displayed.

```
#include <stdio.h>
#include <gsl/gsl_const_mksa.h>

int
main (void)
{
  double c  = GSL_CONST_MKSA_SPEED_OF_LIGHT;
  double au = GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  double minutes = GSL_CONST_MKSA_MINUTE;

  /* distance stored in meters */
  double r_earth = 1.00 * au;
  double r_mars  = 1.52 * au;

  double t_min, t_max;

  t_min = (r_mars - r_earth) / c;
  t_max = (r_mars + r_earth) / c;

  printf ("light travel time from Earth to Mars:\n");
  printf ("minimum = %.1f minutes\n", t_min / minutes);
  printf ("maximum = %.1f minutes\n", t_max / minutes);

  return 0;
}
```

Here is the output from the program,

```
light travel time from Earth to Mars:
minimum = 4.3 minutes
maximum = 21.0 minutes
```

## References and Further Reading

The authoritative sources for physical constants are the 2006 CODATA recommended values, published in the article below. Further information on the values of physical constants is also available from the NIST website.

- P.J. Mohr, B.N. Taylor, D.B. Newell, “CODATA Recommended Values of the Fundamental Physical Constants: 2006”, Reviews of Modern Physics, 80(2), pp. 633–730 (2008).
- <http://www.physics.nist.gov/cuu/Constants/index.html>
- <http://physics.nist.gov/Pubs/SP811/appenB9.html>