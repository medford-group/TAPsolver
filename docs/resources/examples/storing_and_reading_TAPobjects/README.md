
# Storing and Reading TAPobjects

In the previous example, we created a TAPsolver TAPobject for a carbon monoxide oxidation reaction mechanism (named TAP_test). Now, let's brefiely show how we can store and use this information for later.

To store the TAPobject, simply use:

	save_object(TAP_test,'./TAP_test.json')

where the TAPobject TAP_test will be saved as a json file for later use. When desired, the object can be reloaded with:

	read_TAPobject('./TAP_test.json')

This will be very useful when running simulations in parallel or including in publications for submission.