
# Storing and Reading TAPobjects and Synthetic/Experimental Data

In the previous example, we created a TAPsolver TAPobject for a carbon monoxide oxidation reaction mechanism (named TAP_test). Now, let's brefiely show how we can store and use this information for later.

To store the TAPobject, simply use:

	save_object(TAP_test,'./TAP_test.json')

where the TAPobject TAP_test will be saved as a json file for later use. When desired, the object can be reloaded with:

	new_TAP_object = read_TAPobject('./TAP_test.json')

This will be very useful when running simulations in parallel or including in publications for submission.

Last the synthetic data can be read using the TAPsolver function:

	new_data_object = read_experimental_data_object(file_name)

The resulting data object is a dictionary of all the gas fluxes and pulses. For example, if you were pulsing carbon monoxide for 10 pulses and wanted to look at the outlet flux of carbon monoxide over the 7th pulse, you would reference the new_data_object with:

	new_data_object['CO'][7]