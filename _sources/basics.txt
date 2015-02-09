.. _basics:


.. default-domain:: py

Legendre Functions
=======================

.. module:: legendre
   :synopsis: legendre polynomial related functions


.. method:: PlON(lmax, z)

	This function evalutates all ortho-normalized legendre polynomials up to degree ``lmax``.

   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
   :return: A vector of all associated Legendgre polynomials evaluated at ``z`` up to ``lmax``
   :rtype: numpy 1D array


.. method:: PlmON(lmax, z, csphase=1)

	This function evalutates all the normalized associated Legendre
	functions up to degree ``lmax``.

	**Comments from the SHTOOLS library:** The functions are initially scaled by
	10^280 sin^m in order to minimize the effects of underflow at large m
	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
	On a Mac OSX system with a maximum allowable double precision value of
	2.225073858507203E-308 the scaled portion of the algorithm will not overflow
	for degrees less than or equal to 2800.

	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta),
	with the intial value of rescalem being equal to 1/scalef (which is here equal
	to 10^280). This will gradually reduce this huge number to a tiny number, and will
	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
	number will be multipled by the old value of rescalem.

	Temporary variables in saved in an allocated array. In order to explicitly deallocate this
	memory, call this routine with a spherical harmonic degree of -1.


   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
   :param int csphase: Condon-Shortley phase factor

	1  = Do not Apply
	-1 = Apply the phase factor

.. method:: PlmBar(lmax, z, csphase=1)

	This function evalutates all of the normalized associated Legendre
	functions up to degree lmax. 

	**Comments from the SHTOOLS library:** The functions are initially scaled by
	10^280 sin^m in order to minimize the effects of underflow at large m
	near the poles (see Holmes and Featherstone 2002, J. Geodesy, 76, 279-299).
	On a Mac OSX system with a maximum allowable double precision value of
	2.225073858507203E-308 the scaled portion of the algorithm will not overflow
	for degrees less than or equal to 2800.

	For each value of m, the rescaling factor is computed as rescalem=rescalem*sin(theta),
	with the intial value of rescalem being equal to 1/scalef (which is here equal
	to 10^280). This will gradually reduce this huge number to a tiny number, and will
	ultimately underflow. In order to prevent this underflow, when rescalem becomes less than
	10^(-280), the subsequent rescaling factors of sin(theta) will be directly applied to Plm, and then this
	number will be multipled by the old value of rescalem.


   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
   :param int csphase: Condon-Shortley phase factor (1=Do not Apply, -1=Apply)



.. method:: PLegendreA(lmax, z, csphase=1)

   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
   :param int csphase: Condon-Shortley phase factor (1=Do not Apply, -1=Apply)



.. method:: PLegendre(lmax, z)

   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``



.. method:: PlmSchmidt(lmax, z, csphase=1)

   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
   :param int csphase: Condon-Shortley phase factor (1=Do not Apply, -1=Apply)



.. method:: PlSchmidt(lmax, z)

   :param int lmax: Maximum spherical harmonic degree to compute
   :param real z: ``cos(colatitude)`` or ``sin(latitude)``
