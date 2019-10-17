from dessinercestgagner.src.Shape import Shape



if __name__ == "__main__":
    original_shape = Shape('Bone')
    shapes = original_shape.extend()
    for shape in shapes:
        shape.get_profile().display()