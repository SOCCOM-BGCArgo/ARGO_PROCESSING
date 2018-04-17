function y = fresnel(Z)
%  GIVEN ZENITH ANGLE Z, REFLECTANCE R IS CALCULATED USING
%  FRESNEL'S EQUATIONS AS GIVEN BY KIRK (1983).
      NW = 1.341;
      if Z <= 0.0 | Z >= 90.0; error( 'ZENITH ANGLE OUT OF RANGE' ); end
      ACSZ = asin(sin(deg2rad(Z))/NW);
      A = (sin(deg2rad(Z) - ACSZ)).^2 ./ (2.0 *(sin(deg2rad(Z) + ACSZ)).^2);
      B = (tan(deg2rad(Z) - ACSZ)).^2 ./ (2.0 *(tan(deg2rad(Z) + ACSZ)).^2);
      y = A + B;

