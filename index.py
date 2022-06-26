#! python3

import cgitb

cgitb.enable()

# below two lines are compulsory, don't delete them
print("Content-Type: text/html")
print()


# Home page HTML
print(f""" 

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="ie=edge">
    <title>BINF PROJECT</title>
    <link rel="stylesheet" href="./styles.css">
</head>

<body>
   <div class="backdrop"></div>
       <main>
        <form class="input-form" action="./alignSupercompose.py" method="GET">
            <div class="form-control">
                <label for="genename">Gene Name</label>
                <input type="text" name="genename" id="genename" value="MAPT" required>
            </div>
            <div class="form-control">
                <label for="speciesname">Species Name</label>
                <input type="text" name="speciesname" id="speciesname" value="Homo sapiens" required>
            </div>
            
            <div>
                  <input type="checkbox" id="structure" name="structure" value="true">
                  <label for="structure"> Direct Protein Structure</label><br>
            </div>
            <button class="btn" type="submit">Run</button>
        </form>
    </main>
</body>

</html>

 """)
