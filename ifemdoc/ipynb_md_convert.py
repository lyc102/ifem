import os

def find_files(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        for _file in files:
            if name in _file:
                result.append(os.path.join(root, _file))
    return result


if __name__ == '__main__':
    files = find_files('.ipynb', '.')
    for f in files:
        os.system(f'jupyter nbconvert --to markdown {f}')
