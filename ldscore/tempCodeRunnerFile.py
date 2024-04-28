    irwls_fix = IRLS(x, y - irwls.get_intercept(), fix_intercept=True)
    irwls_fix.regression()
    print("h_g^2 fix = %.4f" % irwls_fix.get_coefficients())