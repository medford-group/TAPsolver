class MetaClass(type):    
	def getfoo(self):
		return self._foo
	foo = property(getfoo)

	@property
	def bar(self):
		return self._bar

class MyClass(object, metaclass=MetaClass):
	_foo = 'abc'
	_bar = 'def'

print(MyClass.foo)
print(MyClass.bar)

MyClass._foo = 'test'
print(MyClass.bar)